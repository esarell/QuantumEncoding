'''This program is disigned to run a generic version of the Grover-Rudolph algorithm
It should take a probabiltiy distribution and encoded it into a quantum state of m qubits'''
import qiskit as qt
import matplotlib.pyplot as plt
import numpy as np
import Quantum_Tools as qtool
from qiskit.circuit.library.standard_gates import RYGate
import Linear_Piecewise as lpw
from GenBounds import gen_bounds
import math

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
    This got replaced by the Linear Piece wise function
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
            print(theta_array)
            print(theta_array[m][0])
            theta=theta_array[m][0]
            gate = RYGate(2*theta).control(m)
            qc.append(gate,[*qr[n-m-1:][::-1]])
            for counter,i in enumerate(pattern):
                theta=theta_array[m][counter+1]
                qc.x(qr[n-i:n])
                gate = RYGate(2*theta).control(m)
                qc.append(gate,[*qr[n-m-1:][::-1]])
        return(pattern)    

def Grover_Rudolph_func_small(n,distribution):
    '''Args: This function is for when m<4
    n: the number of qubits for our circuit
    distribution: is the probability distrbution we are trying to encode
    circ:
    qr:
    Returns: an encoded quantum circuit '''
    #Define how many levels we want to have
    m = list(range(0, n-1))
    #m=[0,1]
    #Initalise the quantum circuit
    # our probability distribution will be put on n qubits
    qr= qt.QuantumRegister(n,'q')
    # We then will add 6 bits for the ancillary register 
    #anc = qt.QuantumRegister(6,'ancilla')
    circ = qt.QuantumCircuit(qr)
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
    print(theta_array)
    xGateAdd(len(m)-1,x_gate_pattern,circ,qr,theta_array,n)
    #Draws the quantum circuit
    #qc.draw("mpl")
    #plt.show()
    circ.draw("mpl")
    plt.show()


def ThetaRotation(circ,qr,anc,bit,wrap =True):
    '''Args:
    circ: quantum circuit
    qr: quantum register
    anc: quantum register
    bit: number of which bit in qr register you want rotations applied to 
    wrap: boolean turns indivdual gates into one big gate defualt True 
    Returns: circuit with added Y Rotation gates which are condintional on anc Reg '''
    size_qr = len(qr)
    size_anc =len(anc)
    if wrap == True:
        qr =qt.QuantumRegister(size_qr,"qr")
        anc = qt.QuantumRegister(size_anc,"anc")
        circ = qt.QuantumCircuit(qr,anc)

    for i in range(size_anc):
        temp = i-(size_anc-1)
        theta = 2**(temp)
        rotation_gate = RYGate(theta).control(1)
        circ.append(rotation_gate,[anc[i],qr[bit]])

    if wrap ==True:
        circ = circ.to_gate()
        circ.label = "Ry Theta" 

    return circ
    
def GR_Const_Theta(n,k,theta_array):
    """This is the code which first does general state preparation
    Then after a point k, we do fixed theta values
    Args: n: number of qubits
    k: when you flip systems
    theta: a list of theta values
    Returns: an encoded quantum circuit"""
    fixed_theta=[0.12,0.23,0.4,0.1,0.32,0.64,0.11,0.22]
    control=["000","001","010","011","100","101","110","111"]
    qr= qt.QuantumRegister(n,'q')
    lab=qt.QuantumRegister(1,"L")
    # We then will add 6 bits for the ancillary register 
    anc = qt.QuantumRegister(n,'ancilla')
    circ = qt.QuantumCircuit(qr,lab,anc)
    #Loop through each level of m
    x_gate_pattern =[]    
    xGateAdd(k,x_gate_pattern,circ,qr,theta_array,n)

    amount = n-k-1
    print("amount:",amount)
    #rotation_gate = RYGate(0.3).control(3,ctrl_state="110")
    #circ.append(rotation_gate,[*qr[6:9],qr[4]])
    
    #Adds the fixed theta rotation this if for 3 controls
    #So this is for 8 theta values per each level
    for i in range(amount):
       for n,value in enumerate(fixed_theta):
            rotation_gate = RYGate(value).control(3,ctrl_state=control[n])
            circ.append(rotation_gate,[*qr[6:9],qr[amount-1-i]])

    circ.draw("mpl")
    plt.show()


def GR_function(n):
    """Args: For when m>5
    n,data
    n: the number of qubits for our circuit
    data: 
    We are gonna change data to be a 2d array of all the points from the spile that we want to encode
    Instead of it being a function 
    Returns: an encoded quantum circuit """
    qr= qt.QuantumRegister(size=4,name='q')
    anc = qt.QuantumRegister(size=9,name='anc')
    lab = qt.QuantumRegister(size=3,name='lab')
    target = qt.QuantumRegister(size=1,name='tar')
    coff = qt.QuantumRegister(size=9,name='coff')
    cla_reg =qt.ClassicalRegister(size=9,name="cla")
    circ = qt.QuantumCircuit(qr,target,anc,lab,coff,cla_reg)

    '''Currently this GR will only work for 1 qubit we still need to iterate though it
    Figure out how this works in terms of where the rotations occur, am i looping though it n amount of times in qr
     or for 2**n times. e.g at 3 qubits is it 3 or 8 '''
    #Currently replacing the GR Small function at some point we will only apply this to 
    #bits after m<5 or whatever we decide to do
    initalise_gate = lpw.initalSuperpostion(circ,qr)
    circ.append(initalise_gate,[*qr])
    
    circ.decompose().draw("mpl")
    plt.show()

    '''This will get changed once I am handelling real data '''
    data = gen_bounds(4)
    print(len(data))
    Xdata =data[0]
    Ydata=data[1]
    #Xdata =[1,2,3,4,5,6,7,8,9]
    #Ydata=[3,5,6,7,8,9,9,10,11]
    data = gen_bounds(4)
    print(len(data[0]))
    
    index=[]
    for i in range(4):
        points = (2**i)+1
        #Calculate which datapoints we will need
        if i==0:
            index.append(0)
            index.append(len(data[0])-1)
        else:
            temp =[]
            for place,point in enumerate(index):
                print(point)
                try:
                    half =  int((point+index[place+1])/2)
                    temp.append(point)
                    temp.append(half)
                except:
                    temp.append(point)
            index=temp
        print(index)

        print(points)
    print(data[0][16])
    lpw_gate = lpw.LinearPiecewise(circ,qr,anc,coff,lab,target,Xdata,Ydata)
    circ.append(lpw_gate,[*qr[0],*anc,*coff,*lab,*target])
    
    rotational = ThetaRotation(circ,qr,anc,1)
    circ.append(rotational,[*qr,*anc])


    circ.measure(anc,cla_reg)
    shots = 100
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    counts = result.get_counts(tqc)
    print("counts:",counts)

    circ.decompose().draw("mpl")
    plt.show()
    

    return circ

if __name__ == "__main__":
    #test = qtool.my_binary_repr(1.25,6,nint=1 )
    #print(test)
    #Grover_Rudolph_func_small(4,[0,1,2,3,4,5,6,7,8,10])
    #qr= qt.QuantumRegister(size=4,name='q')
    # We then will add 6 bits for the ancillary register 
    #anc = qt.QuantumRegister(size=9,name='anc')
    #lab = qt.QuantumRegister(size=3,name='lab')
    #target = qt.QuantumRegister(size=1,name='tar')
    #classical = qt.ClassicalRegister(size=4,name="cla")
    #circ = qt.QuantumCircuit(qr,anc)x
    #GR_function(4)
    GR_Const_Theta(9,5,[[0.01],[0.1,0.2],[0.3,0.4,0.5,0.6],[0.07,0.8,0.9,0.10,0.11,0.12,0.13,0.14],[0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31],[0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63]])
    '''
    result = lpw.inputValue(circ,qr,[0,1,1,1]).to_instruction()

    circ.append(result,qr)

    label_Gate_add = lpw.labelGate(circ,qr,target,anc,lab).to_instruction()
    circ.append(label_Gate_add,[*qr,target[0],*anc,*lab,])
    circ.measure(target[0],classical[0])
    shots = 10
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    counts = result.get_counts(tqc)
    print("counts:",counts)'''

    
    #result = inputValue(circ,qr,[1,0,0,0,0,0])

    #lpw.calculateCoffs([4,5],[8,6])


