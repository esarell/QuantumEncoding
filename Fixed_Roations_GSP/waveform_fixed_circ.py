import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
from qiskit.providers.aer import QasmSimulator
import matplotlib.pyplot as plt
import numpy as np
import pickle

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
    
    for count,i in enumerate(condition):
        start = len(i)
        print("i",len(i))
        rotation_gate = RYGate(theta[count]).control(len(i),ctrl_state=i)
        #[*qr[num_qubits-i-1:][::-1]
        #The 8 here is currently hard coded for n=9 qubits, this could be easily changed by passing through n as a parameter
        circ.append(rotation_gate,[*qr[(8-len(i)+1):],qr[qubit_no]])
        #circ.append(rotation_gate,[*qr[9-3:][::-1]])


    if wrap ==True:
        circ = circ.to_gate()
        circ.label = "Ry Theta" 

    return circ

def prob_normalised(data):
    result=[]
    for x in data:
        square = x*x
        norm = np.sqrt(square/sum(data*data))
        result.append(norm)
    return result

def amp_normalised(data):
    result=[]
    for x in data:
        norm = np.sqrt(x/sum(data))
        result.append(norm)
    return result

def waveform_amps():
    with open("waveform_amps.pkl", 'rb') as pickleFile:
        test = pickle.load(pickleFile)
    print(test[0])
    return test

def Fidelity(expected_amps,measured_amps):
    '''
    Args: 
    expected_amps: list of the amplitudes from the discritsed function
    measured_amps: list of amplitudes taken from the statevector of the quantum circuit
    Returns:
    the fidelity (float)
    '''
    current =0
    for count,i in enumerate(expected_amps):
        current = current+(i*measured_amps[count])
    
    fidelity = current*current
    return fidelity


def waveform_theta(frequency,m):
    #[0.74928181,0.73045836,0.72835018,0.71932391,0.73274099,0.74580992,0.75401846]
    j=0
    #This is specifically for f^-7/6 (so p(i)is f^-7/3)
    #Will need changing for a more complicatated distribution
    #Gets a list of normalised amplitudes for all the frequencies 
    temp_amps= waveform_amps()
    amps=prob_normalised(temp_amps)
    thetas=[]
    gsp_theta=[]
    for i in range(m):
        j=pow(2,i)
        if i<3:
            for x in range(j):
                start=int(x*pow(2,m-i))
                mid = int((x+0.5)*pow(2,m-i))
                end = int((x+1)*pow(2,m-i))
                upper = sum(amps[start:(mid)])
                lesser =sum(amps[start:(end)])
                costheta2 = upper/lesser
                costheta = np.sqrt(costheta2)
                theta = np.arccos(costheta) 
                gsp_theta.append(theta*2)
        elif i>=3:
            #print("m:",i)
            if i==3:
                thetas.append(gsp_theta)
                '''Starting points frequency: 40,78.75,117,195,272.5
                #These indexs are hard coded for when n=9 would need to be changed for different points
                #Also these directly correspond with the set up of the circuit so would also need to be changed if the circuit changes'''
                index=[0,64,128,256,384]
                #index=[0,4,8,12,16,24]
            else:
                '''Set starting points for all the
                40, 59.375, 78.75, 117.5, 195, 272.5
                0,32,64,128,256,384'''
                index=[0,32,64,128,256,384]
                #index=[0,2,4,8,16,24]
            temp_theta =[]
            for current_index in index:
                temp_start=current_index
                print("start:",frequency[temp_start])
                place = current_index/pow(2,m-i)
                print("x:",place)
                increment = int(pow(2,m)/pow(2,i))
                mid =int((place+0.5)*pow(2,m-i))
                end=int((place+1)*pow(2,m-i))
                print("increment:",increment)
                print("mid",mid)
                print("end:",end)
                upper = sum(amps[current_index:mid])
                lesser =sum(amps[current_index:end])
                costheta2 = upper/lesser
                costheta = np.sqrt(costheta2)
                theta = np.arccos(costheta) 
                print("theta:",theta)
                temp_theta.append(theta*2)
            thetas.append(temp_theta)
    if i<3:
        thetas.append(gsp_theta)
    #print("theata:",thetas)
    return thetas


def Waveform_Fixed_Rots(n):
    """Ashwins theta values
    [0.69134909] 
 [0.75468131 0.498768  ] 
 [0.7304692  0.75837712 0.70859427 0.66545498] 
 [0.72251893 0.74890938 0.75837712 0.70859427 0.66545498] 
 [0.73581314 0.74934639 0.74890938 0.75837712 0.70859427 0.66545498] """
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)
    #basic = General_State_Prep([0.69134909,0.75468131,0.498768,0.7304692,0.75837712, 0.70859427, 0.66545498])
        #Generate theta values for the GSP section
    #A list of frequency values evenly spaced, the amount is based on n
    frequency = np.linspace(40,350,num=pow(2,n),endpoint=False)
    angles = waveform_theta(frequency,n)
    basic = General_State_Prep(angles[0])
    
    circ.append(basic,qr[n-3:n])
    circ.save_statevector(label='v1')
    #This is for m=3
    #circ.save_statevector()
    controls =["000","001","01","10","11"]
    m_3=ThetaRotation(circ,qr,controls,5,angles[1],True)
    circ.append(m_3,[*qr])
    #for the rest of the levels m>3
    circ.save_statevector(label='v2')
    
    for i in range((n-4)):
        controls =["0000","0001","001","01","10","11"]
        #theta_m4=[0.73581314,0.74934639,0.74890938,0.75837712,0.70859427,0.66545498]
        #theta_m4=[[0.74066914,0.75075189,0.75726406,0.76502399,0.77230382,0.77575838],[0.75760343,0.76520445,0.76942786,0.77434713,0.77851473,0.78040074],[0.76978013,0.77445391,0.77707962,0.77963685,0.78186772,0.78285331],[0.77710115,0.77969497,0.78110797,0.78245585,0.78361015,0.78411397],[0.78111943,0.78248617,0.78321923,0.77277432,0.78449838,0.7847531]]
        fixed=ThetaRotation(circ,qr,controls,8-(4+i),angles[i+2],True)
        circ.append(fixed,[*qr])
    circ.save_statevector(label='v3')
    backend = QasmSimulator()
    backend_options = 'statevector'
    job = qt.execute(circ, backend, shots=1000)
    result = job.result()
    statevector=result.data(0)['v1']
    statevector2=result.data(0)['v2']
    statevector3=result.data(0)['v3']
    circ.decompose().draw("mpl",fold=-1)
    plt.show()
    print(statevector)
    print(statevector2)
    print(statevector3)

    #amps = waveform_amps()
    temp_amps= waveform_amps()
    amps=amp_normalised(temp_amps)
    #amps_7_6 = amplitudes_7_over_6(frequency)
    plt.plot(frequency,np.sqrt(statevector3.probabilities()),color='k',label='Statevector')
    plt.plot(frequency,amps,color='r',label='Waveform')
    #plt.plot(frequency,amps_7_6,color='b',label='Amps -7/6')
    plt.legend()
    plt.xlabel('f (Hz)')
    plt.ylabel('A')
    plt.show()
    print("Fidelity: ",Fidelity(amps,np.sqrt(statevector3.probabilities())))



if __name__ == "__main__":
    waveform_amps()
    Waveform_Fixed_Rots(9)