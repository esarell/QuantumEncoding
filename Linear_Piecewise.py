import qiskit as qt
import matplotlib.pyplot as plt
import numpy as np

def calculateCoffs(X_data,Y_data):
    """
    Args:
    Takes in two data points
    Returns:
    List with two coffcents 
    """
    A = np.vstack([X_data, np.ones(len(X_data))]).T
    print(A)
    m = (Y_data[0]-Y_data[1])/(X_data[0]-X_data[1])
    c = -(m*X_data[0]) + Y_data[0]

    result =[m,c]
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
    #Define ncut

    #You can change step size if we want smaller bounds
    for i,bound in enumerate(range(0,bound)):
        
        intcomp_gate = qtool.integer_compare(circ, qr, target, anc, bound, wrap=True, uncomp=False, label='P'+str(i))
        circ.append(intcomp_gate, [*qr, target[0], *anc[:]])

        inc_gate = qtool.increment_gate(circ, lab, wrap=True, label='SET'+str(i), ncut=0, QFT_on=False, iQFT_on=False).control(1)
        circ.append(inc_gate, [target[0], *lab[:]])

        intcomp_gate_inv = qtool.integer_compare(circ, qr, target, anc, bound, wrap=True, uncomp=False, inverse=True, label='P'+str(i))
        circ.append(intcomp_gate_inv, [*qr, target[0], *anc[:]])
        
    QFT_Gate = qtool.QFT(circ,lab,wrap=True)
    circ.append(QFT_Gate,lab)  
    #Left over from fergus code unsure why this here currently
    #circ.x(qr[-1])
    return circ

def LinearPiecewise(circ,lab,coeffiecents):
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
        #NOTE: These are currently hardcoded this needs to be generalised
        #Easy to do, just need to make the number be variables based on the size of register being sent in
        qr = qt.QuantumRegister(6, 'q_reg')
        target = qt.QuantumRegister(1, 'q_targ')
        coff = qt.QuantumRegister(6,"q_coff")
        anc = qt.QuantumRegister(6, 'q_ans')
        lab = qt.QuantumRegister(5, 'q_lab')
        circ = qt.QuantumCircuit(qr, target, anc, lab) 

    input_gate_add =inputValue(circ,qr,[1,0,0,0,0,0])
    circ.append(input_gate_add,qr)
    #Adds the labels to the different subdomains based on the value from 
    label_Gate_add = labelGate(circ,qr,anc,lab,target)
    circ.append(label_Gate_add,[*qr,*anc,*lab,target[0]])
    #Load in coefficents
    #Convert coeffiencents to binary representation


    input_gate_add =inputValue(circ,coff,[1,0,0,0,0,0])

    #Multiply into the x register

    #Unload coefficents

    #Load in coefficents

    #Add the coefficent with 

    #Unload coefficents


    #circ.draw("mpl")
    #plt.show()

if __name__ == "__main__":
    #test = qtool.my_binary_repr(1.25,6,nint=1 )
    #print(test)
    '''qr= qt.QuantumRegister(size=4,name='q')
    # We then will add 6 bits for the ancillary register 
    anc = qt.QuantumRegister(size=4,name='anc')
    lab = qt.QuantumRegister(size=3,name='lab')
    target = qt.QuantumRegister(size=1,name='tar')
    classical = qt.ClassicalRegister(size=4,name="cla")
    circ = qt.QuantumCircuit(qr,anc,lab,target,classical)
    result = inputValue(circ,qr,[1,0,0,0])
    circ.append(result,qr)
    label_Gate_add = labelGate(circ,qr,anc,lab,target)
    circ.append(label_Gate_add,[*qr,*anc,*lab,target[0]])
    #result = inputValue(circ,qr,[1,0,0,0,0,0])
    circ.draw("mpl")
    plt.show()'''
    calculateCoffs([4,5],[8,6])