import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
import matplotlib.pyplot as plt

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
    #xGateAdd(k,x_gate_pattern,circ,qr,theta_array,n)

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

#This is Roselyns code
def GenBinaryStrings(length):
    binary_strings = []
    max_num = 2 ** length

    for i in range(max_num):
        binary_string = format(i, 'b').zfill(length)
        binary_strings.append(binary_string)
    return binary_strings


#This is Roselyns code
def ReverseString(stringlist):
    reversed_list = [string[::-1] for string in stringlist]
    return reversed_list

#This is Roselyns code with my additions
def General_State_Prep(rotations):
    index=1
    num_qubits=len(rotations).bit_length()
    #builds blank circuit
    circuit = qt.QuantumCircuit(num_qubits)
    #first rotation is the first item in rotations so can just be directly placed
    circuit.ry(rotations[0],0)

    #generates binaries
    for i in range(1,num_qubits):
        binaries = GenBinaryStrings(i)
        binaries = ReverseString(binaries)
        qbits = list(range(0,i+1))
        #use binaries to construct the controlled y rotations
        for string in binaries:
            control_y = RYGate(rotations[index]).control(i,None,string)
            circuit.append(control_y,qbits)
            index += 1

    return circuit



def ThetaRotation(circ,qr,condition,qubit_no,theta,wrap =True):
    '''Args:
    circ: quantum circuit
    qr: quantum register
    condition:list of str what conditions should it have
    qubit_no:which qubit
    wrap: boolean turns indivdual gates into one big gate defualt True 
    Returns: '''
    size_qr = len(qr)
    if wrap == True:
        qr =qt.QuantumRegister(size_qr,"qr")
        circ = qt.QuantumCircuit(qr)
    
    for i in condition:
        start = len(i)
        rotation_gate = RYGate(theta).control(len(i),ctrl_state=i)
        circ.append(rotation_gate,[*qr[0:len(i)],qr[qubit_no]])


    if wrap ==True:
        circ = circ.to_gate()
        circ.label = "Ry Theta" 

    return circ


def Inspiral_Fixed_Rots(n):
    """Create a circuit for the inspiral f^-7/6
    Args: n: number of qubits
    Returns: Circuit """
    qr= qt.QuantumRegister(size=n,name='q')
    circ = qt.QuantumCircuit(qr)
    basic = General_State_Prep([0.1,0.2,0.3,0.4,0.5,0.6,0.7])
    circ.append(basic,qr[0:3])
    #This is for m=3
    controls =["000","100","010","110","01","11"]
    m_3=ThetaRotation(circ,qr,controls,3,0.1,True)
    circ.append(m_3,[*qr])
    #for the rest of the levels m>3
    for i in range((n-4)):
        controls =["0000","1000","100","10","01","11"]
        fixed=ThetaRotation(circ,qr,controls,(4+i),0.1,True)
        circ.append(fixed,[*qr])
    circ.decompose().draw("mpl",fold=-1)
    plt.show()

def Waveform_Fixed_Rots(n):
    """ """
    qr= qt.QuantumRegister(size=n,name='q')
    circ = qt.QuantumCircuit(qr)
    basic = General_State_Prep([0.1,0.2,0.3,0.4,0.5,0.6,0.7])
    circ.append(basic,qr[0:3])
    #This is for m=3
    controls =["000","100","10","01","11"]
    m_3=ThetaRotation(circ,qr,controls,3,0.1,True)
    circ.append(m_3,[*qr])
    #for the rest of the levels m>3
    for i in range((n-4)):
        controls =["0000","1000","100","10","01","11"]
        fixed=ThetaRotation(circ,qr,controls,(4+i),0.1,True)
        circ.append(fixed,[*qr])
    circ.decompose().draw("mpl",fold=-1)
    plt.show()


if __name__ == "__main__":
    
    Inspiral_Fixed_Rots(9)
    Waveform_Fixed_Rots(9)