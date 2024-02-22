'''This program is disigned to run a generic version of the Grover-Rudolph algorithm
It should take a probabiltiy distribution and encoded it into a quantum state of m qubits'''
import qiskit as qt
import matplotlib.pyplot as plt
import numpy as np
import scipy
import Quantum_Tools as qtool
from qiskit.circuit.library.standard_gates import RYGate

def prob(qubit,fmin,fdelta,distrbution_type):
    '''
    Args:xs
    qubit: takes in the Qubit (we will use this in a decimal form)
    fmin: minimum frequency we are using
    fdelta: the change in frequence per each step
    distribution_type: the number of the power to be applied to this
    Note in Hayes paper there distribution_type is 7/3
    Return:
    probabilty from that point in the qubit :p(i) in the 
    '''
    '''f=fmin + i delta f
    so p(i) = (fmin+i delta f)^(=7/3) / Nc
    Where Nc is the normalisation constant
    delta f = (fmin-fmax)/2^n '''
    #print(qubit)
    probabilityI = (fmin + (qubit*fdelta))**distrbution_type
    #This is where we need the probability distribution
    #print("prob",probabilityI)
    return probabilityI

def cos2_theta(m,j,n,distribution):
    '''Args:
    m: the level we are at
    j: which bin we are in
    n:
    distribution:
    Return:
    cos^2 theat(m,j)
    '''
    #two sums which are then divided over each other
    i = j*2**(n-m-1)
    #I think this is where i need to minus 1 from the equation
    #Think about it for the sum though?
    half = (j+0.5)*2**(n-m-1)
    full =  (j+1)*2**(n-m-1)
    #print(i)
    #Calculate the constants from the classical distribiution for prob
    fmin = distribution[0]
    size = len(distribution)
    fmax = distribution[(size-1)]
    #normalise our distribution
    normal = np.sqrt(np.sum(np.abs(distribution)**2))
    fdelta = (fmax-fmin)/(2**n)
    top =0
    bottom =0
    #CHANGE FOR THE DIFFERENT KINDS OF PROBABILITY
    #THIS CURRENTLY DOES NOTHING
    distrbution_type =1
    for x in range(i,int(half)):
        top = top + prob(x,fmin,fdelta,distrbution_type)
    for y in range(i,full):
        bottom = bottom + prob(y,fmin,fdelta,distrbution_type)
    result = top/bottom
    #Now get the angle from the cos^2(theta)
    result = np.arccos(np.sqrt(result))
    return result
def Q_operator(m,j,theta,circ,qr,n):
    '''
    Args:
    m = level
    j = bin
    theta= outcome from the cos2_theta fucntion
    circ= current quantum circuit
    Return:
    circ = Updated quantum circuit with the changed ancillary register
    '''
    #First convert from decimal to binary using qiskit tools
    test = qtool.my_binary_repr(theta,6,nint=1 )
    if m == 0:
        #level 1 only requires one rotational gate
        circ.ry(2*theta,qr[0:n-1])
    return test

def setAncillary(theta,qubit,circ,anc):
    """
    theta:
    qubit: the number we are currently on in the loop this can be used for our control bit?
    circ:
    anc:
    return

    """

    thetaBinary = qtool.my_binary_repr(theta,6,nint=1 )
    #print('anc',len(anc))
    for i in range(len(anc)):
        if thetaBinary[i] =='1':
            circ.x(anc[i])

def inputValue(circ,qr,value):
    """Args:
    circ: Current circuit
    qr: The register that we want to put our value in
    value: the number we want to encode (in binary)
    Return:
    Circuit that has the value in the qr 
    This is just computation bais encoding"""
    if True:
        qr = qt.QuantumRegister(4, 'q_reg')
        circ = qt.QuantumCircuit(qr) 
    
    for i,value in enumerate((value[::-1])):
        print(i)
        if value == 1:
            circ.x(qr[i])
            print("appended at qubit:",i)
    #circ.draw("mpl")
    #plt.show()
    return circ

def labelGate(circ,qr,anc,lab,target):
    """
    Args:
    circ:
    qr: register containing the value x
    anc:
    lab:
    target:

    Return:
    Goal is to put out labels into the lab register
    """

    #We want to Wrap the stuff we are doing gets werid when we don't
    if True:
        qr = qt.QuantumRegister(4, 'q_reg')
        target = qt.QuantumRegister(1, 'q_targ')
        anc = qt.QuantumRegister(4, 'q_ans')
        lab = qt.QuantumRegister(3, 'q_lab')
        circ = qt.QuantumCircuit(qr, target, anc, lab) 
    n=len(qr)
    #Left over from fergus code unsure why this here currently
    #circ.x(qr[-1])
    #Adding QFT to do counting in phase space
    QFT_Gate = qtool.QFT(circ,lab,wrap=True)
    circ.append(QFT_Gate,lab)
    #We should define bounds, currently pic a constant
    bound =2**(n-1)
    print("bound:",bound)
    #Define ncut

    #You can change step size if we want smaller bounds
    for i,bound in enumerate(range(0,bound)):
        
        intcomp_gate = qtool.integer_compare(circ, qr, target, anc, bound, wrap=True, uncomp=False, label='P'+str(i))
        circ.append(intcomp_gate, [*qr, target[0], *anc[:]])

        inc_gate = qtool.increment_gate(circ, lab, wrap=True, label='SET'+str(i), ncut=0, QFT_on=False, iQFT_on=False).control(1)
        circ.append(inc_gate, [target[0], *lab[:]])

        intcomp_gate_inv = qtool.integer_compare(circ, qr, target, anc, bound, wrap=True, uncomp=False, inverse=True, label='P'+str(i))
        circ.append(intcomp_gate_inv, [*qr, target[0], *anc[:]])

        print("bound",bound)
        print("i:",i)
        
    QFT_Gate = qtool.QFT(circ,lab,wrap=True)
    circ.append(QFT_Gate,lab)  
    #circ.x(qr[-1])

    return circ

def LinearPiecewise(circ,lab):
    """
    Function 
    Args:
    circ: the overall circit
        *Note circ will be made of 4 registers that are used and 1 register that isn't used by this function
        
        anc:
        lab: Label Gate
        -------------------------------------------
        qr: This register the main register used by the rest of the grover-rudolph algorithm

    Returns:
    the function f(x) in the output register
    """
    #Split the function into sections calculate the thing it's doing off that classically
    #For this would this be my cos^2 theta 
    #Label register: label which section we are in
    #Coefficent register: Store coeefficent control on label
    if True:
        qr = qt.QuantumRegister(6, 'q_reg')
        target = qt.QuantumRegister(1, 'q_targ')
        anc = qt.QuantumRegister(6, 'q_ans')
        lab = qt.QuantumRegister(5, 'q_lab')
        circ = qt.QuantumCircuit(qr, target, anc, lab) 

    input_gate =inputValue(circ,qr,[1,0,0,0,0,0])
    circ.append(input_gate,qr)
    label_gate_add =labelGate(circ,qr,anc,lab,target)
    circ.append(label_gate_add,lab)
    circ.draw("mpl")
    plt.show()

def xGateAdd(m,pattern,qc,qr,theta_array,n):
    '''Adds the X gates and the Rcy Gates for small values of m
        This is a recursive function
    Args:
    m the max level for the Grover Rudolph Algorithm
    pattern: a blank list which will be populated with the pattern of the X gates
    qc: quantum circuit
    qr: quantum register this is being added to
    theta_array: a 2 dimensional array containing the values for the Rcy gate
    n: the number of qubits for the qr register
    Returns:
    patterns: A list of the X-gate pattern
     '''
    if m == 0:
        theta = theta_array[0][0]
        qc.ry(2*theta,qr[n-1])
        return []
    else:
        x = m -1
        pattern = xGateAdd(x,pattern,qc,qr,theta_array,n)
        temparray = pattern.copy()
        pattern.append(m)
        if m ==1:
            qc.x(n-1)
            theta=theta_array[m][0]
            gate = RYGate(2*theta).control(m)
            qc.append(gate,[n-1,n-2])
            qc.x(n-1)
            theta=theta_array[m][1]
            gate = RYGate(2*theta).control(m)
            qc.append(gate,[n-1,n-2])
            return(pattern)
        else:
            pattern.extend(temparray)
            qc.x(qr[n-m:n])
            theta = theta=theta_array[m][0]
            gate = RYGate(2*theta).control(m)
            qc.append(gate,[*qr[n-m-1:][::-1]])
            for counter,i in enumerate(pattern):
                theta=theta_array[m][counter+1]
                qc.x(qr[n-i:n])
                gate = RYGate(2*theta).control(m)
                qc.append(gate,[*qr[n-m-1:][::-1]])
        return(pattern)    

def Grover_Rudolph_func(n,distribution):
    '''Args: 
    n: the number of qubits for our circuit
    distribution: is the probability distrbution we are trying to encode
    Returns: an encoded quantum circuit '''
    #Define how many levels we want to have
    m = list(range(0, n-1))
    #m=[0,1]
    #Initalise the quantum circuit
    # our probability distribution will be put on n qubits
    qr= qt.QuantumRegister(n,'q')
    # We then will add 6 bits for the ancillary register 
    anc = qt.QuantumRegister(6,'ancilla')
    qc = qt.QuantumCircuit(qr,anc)
    theta_array =[]
    #Loop through each level of m
    x_gate_pattern =[]    
    for i in m:
        if i > 4:
            break
        #split up the probability distribution into two parts
        #This is our j in the maths
        #Define the number of bins j (how many bins we split the probability distribution into)
        current_bins = list(range(0,2**i))
        temp =[]
        for j in current_bins:
            #Calculates the angle based on the probability of that bin m,j
            theta=cos2_theta(i,j,n,distribution)
            temp.append(theta)
            #place=str(i)+str(j)
            #angles[place]= theta
        theta_array.append(temp)
    xGateAdd(len(m)-1,x_gate_pattern,qc,qr,theta_array,n)
    #Draws the quantum circuit
    qc.draw("mpl")
    plt.show()

if __name__ == "__main__":
    #test = qtool.my_binary_repr(1.25,6,nint=1 )
    #print(test)
    #Grover_Rudolph_func(5,[0,1,2,3,4,5,6,7])
    qr= qt.QuantumRegister(size=4,name='q')
    # We then will add 6 bits for the ancillary register 
    anc = qt.QuantumRegister(size=4,name='anc')
    lab = qt.QuantumRegister(size=3,name='lab')
    target = qt.QuantumRegister(size=1,name='tar')
    classical = qt.ClassicalRegister(size=4,name="cla")
    circ = qt.QuantumCircuit(qr,anc,lab,target,classical)
    result = inputValue(circ,qr,[1,0,0,0])
    circ.append(result,qr)
    labelGate(circ,qr,anc,lab,target)
    #result = inputValue(circ,qr,[1,0,0,0,0,0])

    '''circ.append(result,qr)
    circ.measure(qr,classical)
    shots = 10
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    counts = result.get_counts(tqc)
    print("counts:",counts) '''

