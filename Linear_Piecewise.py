import qiskit as qt
import matplotlib.pyplot as plt
import numpy as np
import Quantum_Tools as qtool

def quick_measure(circ,q_reg,cla_reg):
    '''Just a simple function that will do a measurement 10 times and then print the results
    Args:
    circ: overal circuit
    q_reg: the quantum register we are measureing
    cla_reg: the classical register we are putting the result into
    Return: 
    counts: number  '''
    circ.measure(q_reg,cla_reg)
    shots = 10
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    counts = result.get_counts(tqc)
    print("counts:",counts)
    return counts

def calculateCoffs(X_data,Y_data):
    """
    Args:
    Takes in two data points
    Returns:
    List with two coffcents 
    """
    m = (Y_data[0]-Y_data[1])/(X_data[0]-X_data[1])
    c = -(m*X_data[0]) + Y_data[0]

    result =[c,m]
    print(result)
    return result

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
        if value == 1:
            circ.x(qr[i])
            print("appended at qubit:",i)
    #circ.draw("mpl")
    #plt.show()
    return circ

def labelGate(circ,qr,target,anc,lab):
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
        #classical = qt.ClassicalRegister(1,"classical")
        circ = qt.QuantumCircuit(qr, target, anc, lab) 
    n=len(qr)
    #Left over from fergus code unsure why this here currently
    #circ.x(qr[-1])
    #Adding QFT to do counting in phase space
    QFT_Gate = qtool.QFT(circ,lab,wrap=True)
    circ.append(QFT_Gate,lab)
    #We should define bounds, currently pic a constant
    bound =7
    #Define ncut

    #You can change step size if we want smaller bounds
    for i,bound in enumerate(range(0,bound)):
        print("bound",bound)
        #quick_measure(circ,qr,classical)
        intcomp_gate = qtool.integer_compare(circ, qr, target, anc, bound, uncomp=False, label='P'+str(i))
        circ.append(intcomp_gate, [*qr, target[0], *anc[:]])

        inc_gate = qtool.increment_gate(circ, lab, wrap=True, label='SET'+str(i), ncut=0, QFT_on=False, iQFT_on=False).control(1)
        circ.append(inc_gate, [target[0], *lab[:]])

        intcomp_gate_inv = qtool.integer_compare(circ, qr, target, anc,bound, wrap=True, uncomp=False, inverse=True, label='P'+str(i))
        circ.append(intcomp_gate_inv, [*qr, target[0], *anc[:]])
  
    QFT_Gate_inv = qtool.QFT(circ,lab,wrap=True,inverse=True,do_swaps=False)
    circ.append(QFT_Gate_inv,lab)

    #Left over from fergus code unsure why this here currently
    #circ.x(qr[-1])

    if True:
        circ = circ.to_gate()
        circ.label = "LABEL GATE"

    return circ

def load_coefficents(circ, qcoff, qlab, coeffs_in, nint=None, phase=False, wrap=False, inverse=False, label='Load\nCoff', comp2=True):
    """"Adapted from first_gate as seen in the qiskit_tools.py file"""
    n = len(qcoff)
    nlab = len(qlab)
    print("n: ",n)
    print("nlab: ",nlab)
    if True:
        qcoff = qt.QuantumRegister(n, 'q_coff')
        qlab = qt.QuantumRegister(nlab, 'q_lab')
        circ = qt.QuantumCircuit(qcoff, qlab)  

    #Loops though the coefficents
    for i in np.arange(len(coeffs_in)):
        #converts index to binary
        print("label i:",i)
        control_bits = qtool.my_binary_repr(i, nlab, nint=None, phase=False)
        print(control_bits)
        #if greater than one takes the previous number and does in the binary but splitting the highest and lowest order
        if i>0:
            prev_control = qtool.my_binary_repr(i-1, nlab, nint=nlab, phase=False)[::-1]
        else:
            #if eq or less than zero then sets it all to zero
            prev_control = np.ones(nlab).astype(int).astype(str)
        
        #Compares the current control bit with the previous controll bit
        for j,control_bit in enumerate(control_bits[::-1]):
            #Why? flip the bit if its an 1 and the previous was a zero
            #This would increment the value of the lab
            if control_bit=='0' and prev_control[j]=='1':
                circ.x(qlab[j])
        
        if comp2:
            #input_gate = inputValue(, circ, wrap=True).control(nlab)
            current_coff =qtool.my_binary_repr(coeffs_in[i], n, nint=nint, phase=phase)
            input_gate =inputValue(circ, qcoff,current_coff).control(nlab)
            circ=circ.compose(input_gate, [*qcoff,*qlab])
        '''else:
            input_gate = inputValue(qtool.my_binary_repr(np.abs(coeffs_in[i]), n-1, nint=nint, phase=False), circ, reg=qcoff[:-1], wrap=True).control(nlab)
            circ.append(input_gate, [*qlab, *qcoff[:-1]]);
            if coeffs_in[i]<0.:
                circ.append(XGate().control(nlab), [*qlab, qcoff[-1]]);'''

        if i<len(coeffs_in)-1:
            prev_control = qtool.my_binary_repr(i+1, nlab, nint=nlab, phase=False)[::-1]
        else:
            prev_control = np.ones(nlab).astype(int).astype(str)
        
        for j,control_bit in enumerate(control_bits[::-1]):
            if control_bit=='0' and prev_control[j]=='1':
                #Undos the increment of the lab
                circ.x(qlab[j])
    if True:
        circ = circ.to_gate()
        circ.label = label
        if inverse:
            circ = circ.inverse()
            circ.label = label+'†'
        
    return circ

def LinearPiecewise(circ,qr,anc,coff,lab,target,Xdata,Ydata):
    """
    Function 
    Args:
    circ: the overall circit
        *Note circ will be made of 4 registers that are used and 1 register that isn't used by this function
        
        anc:
        lab: lable register
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
        #NOTE: These are currently hardcoded this needs to be generalised
        #Easy to do, just need to make the number be variables based on the size of register being sent in
        qr = qt.QuantumRegister(4, 'q_reg')
        #At some point make target bit be the first bit in the coff register
        #reduce the ampunt og qubits we are using
        anc = qt.QuantumRegister(4, 'q_ans')
        coff = qt.QuantumRegister(4,"q_coff")
        lab = qt.QuantumRegister(3, 'q_lab')
        target = qt.QuantumRegister(1, 'q_targ')
        circ = qt.QuantumCircuit(qr,anc,coff,lab,target) 

    #Bounds are gonna be equal to the Xdata?
    bounds = Xdata
    #Generate coefficents
    A1_coeffs =[]
    A0_coeffs = []
    for i,value in enumerate(Xdata[:-1]):
        print(i)
        try:
            Coeffs =  calculateCoffs(Xdata[i:i+2],Ydata[i:i+2])
        except:
            Coeffs =  calculateCoffs(Xdata[i:],Ydata[i:])
        A0_coeffs.append(Coeffs[0])
        A1_coeffs.append(Coeffs[1])
    
    print("A1",A1_coeffs)
    print("A0",A0_coeffs)

    input_gate_add =inputValue(circ,qr,[1,0,0,0])
    circ.append(input_gate_add,qr)
    #Adds the labels to the different subdomains based on the value from 
    label_Gate_add = labelGate(circ,qr,target,anc,lab)
    circ.append(label_Gate_add,[*qr,target[0],*anc,*lab,])

    #Convert coeffiencents to binary representation
    input_gate_add = load_coefficents(circ,coff,lab,A1_coeffs,wrap=True)
    circ.append(input_gate_add,[*lab,*coff])
    #Multiply with the x register into anc
    multiple_gate_add = qtool.QFTMultiply(circ,qr,coff,anc,nint3=4,wrap=True)
    circ.append(multiple_gate_add,[*qr,*coff,*anc])

    #Unload coefficents
    input_gate_add = load_coefficents(circ,coff,lab,A1_coeffs,inverse=True)
    circ.append(input_gate_add,[*lab,*coff])
    #Load in coefficents

    #Add the coefficent with 

    #Unload coefficents


    #circ.draw("mpl")
    #plt.show()
    return circ

if __name__ == "__main__":
    #test = qtool.my_binary_repr(1.25,6,nint=1 )
    #print(test)
    qr= qt.QuantumRegister(size=4,name='q')
    # We then will add 6 bits for the ancillary register 
    anc = qt.QuantumRegister(size=4,name='anc')
    lab = qt.QuantumRegister(size=3,name='lab')
    target = qt.QuantumRegister(size=1,name='tar')
    coff = qt.QuantumRegister(size=4,name='coff')
    cla_reg =qt.ClassicalRegister(size=3,name="cla")
    circ = qt.QuantumCircuit(qr,target,anc,lab,coff,cla_reg)
    '''qr= qt.QuantumRegister(size=4,name='q')
    # We then will add 6 bits for the ancillary register 
    anc = qt.QuantumRegister(size=4,name='anc')
    lab = qt.QuantumRegister(size=3,name='lab')

    target = qt.QuantumRegister(size=1,name='tar')
    classical = qt.ClassicalRegister(size=4,name="cla")
    circ = qt.QuantumCircuit(qr,anc,coff,lab,target,classical)'''
    Xdata =[1,2,6,9,10,12,13,14,15]
    Ydata=[2,3,6,7,9,12,14,17,18]
    lpw_gate = LinearPiecewise(circ,qr,anc,coff,lab,target,Xdata,Ydata)
    circ = circ.compose(lpw_gate,[*qr,*anc,*coff,*lab,*target])
    #gate_1 = qtool.QFTMultiply(circ,qr,coff,anc,nint3=4,wrap=True)
    #circ = circ.compose(gate_1,[*qr,*coff,*anc])
    '''result = inputValue(circ,qr,[1,0,0,0])
    circ.append(result,qr)
    label_Gate_add = labelGate(circ,qr,anc,lab,target)
    circ.append(label_Gate_add,[*qr,*anc,*lab,target[0]])
    #result = inputValue(circ,qr,[1,0,0,0,0,0])
    '''
    #input_gate_add = load_coefficents(circ,coff,lab,Xdata).to_gate()
    #circ.append(input_gate_add,[*lab,*coff])
    #calculateCoffs([4,5],[8,6])
    circ.measure(lab,cla_reg)
    shots = 100
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    counts = result.get_counts(tqc)
    print("counts:",counts)

    circ.draw("mpl")
    plt.show()