import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
from qiskit.providers.aer import QasmSimulator
import matplotlib.pyplot as plt
import numpy as np

def GenBinaryStrings(length):
    """

    :param length:
    :return:
    """
    binary_strings = []
    max_num = 2 ** length

    for i in range(max_num):
        binary_string = format(i, 'b').zfill(length)
        binary_strings.append(binary_string)
    return binary_strings


#This is Roselyns code
def ReverseString(stringlist):
    """

    :param stringlist:
    :return:
    """
    reversed_list = [string[::-1] for string in stringlist]
    return reversed_list

#This is Roselyns code with my additions
def General_State_Prep(rotations):
    """

    :param rotations:
    :return:
    """
    print("rotations:")
    print(rotations)
    index=1
    num_qubits=len(rotations).bit_length()
    #builds blank circuit
    qr =qt.QuantumRegister(num_qubits,"qr")
    circuit = qt.QuantumCircuit(qr)
    #first rotation is the first item in rotations so can just be directly placed
    circuit.ry(rotations[0],num_qubits-1)

    #generates binaries
    for i in range(1,num_qubits):
        binaries = GenBinaryStrings(i)
        binaries = ReverseString(binaries)
        qbits = list(range(0,i+1))
        #use binaries to construct the controlled y rotations
        for string in binaries:
            control_y = RYGate(rotations[index]).control(i,None,string)
            circuit.append(control_y,[*qr[num_qubits-i-1:][::-1]])
            index += 1

    return circuit


def ThetaRotation(circ,qr,condition,qubit_no,theta_values,wrap =True):
    """

    :param circ: quantum circuit
    :param qr: quantum register
    :param condition: list of str what conditions should it have
    :param qubit_no:
    :param theta:
    :param wrap: boolean turns indivdual gates into one big gate defualt True
    :return:
    """

    size_qr = len(qr)
    print("size_qr",size_qr)
    if wrap == True:
        qr =qt.QuantumRegister(size_qr,"qr")
        circ = qt.QuantumCircuit(qr)
    print(len(condition))
    print("theats",len(theta_values))
    for count,i in enumerate(condition):
        start = len(i)
        rotation_gate = RYGate(theta_values[count]).control(len(i),ctrl_state=i)
        #[*qr[num_qubits-i-1:][::-1]
        #The 8 here is currently hard coded for n=8 qubits, this could be easily changed by passing through n as a parameter
        print("qubit:",(4-len(i)+1))
        circ.append(rotation_gate,[*qr[(4-len(i)+1):],qr[qubit_no]])
        #circ.append(rotation_gate,[*qr[9-3:][::-1]])


    if wrap ==True:
        circ = circ.to_gate()
        circ.label = "Ry Theta"

    return circ

def prob_normalised(data):
    """

    :param data:
    :return:
    """
    data = np.array(data)
    result=[]
    for x in data:
        square = x*x
        norm = np.sqrt(square/sum(data*data))
        result.append(norm)
    return result

def amp_normalised(data):
    """

    :param data:
    :return:
    """
    result=[]
    for x in data:
        norm = np.sqrt(x/sum(data))
        result.append(norm)
    return result

def Fidelity(expected_amps,measured_amps):
    """

    :param expected_amps: list of the amplitudes from the discritsed function
    :param measured_amps: list of amplitudes taken from the statevector of the quantum circuit
    :return: the fidelity (float)
    """
    current =0
    for count,i in enumerate(expected_amps):
        current = current+(i*measured_amps[count])

    fidelity = current*current
    return fidelity