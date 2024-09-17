import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
from qiskit.providers.aer import QasmSimulator
import matplotlib.pyplot as plt
import csv
import numpy as np

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
        rotation_gate = RYGate(theta[count]).control(len(i),ctrl_state=i)
        circ.append(rotation_gate,[*qr[0:len(i)],qr[qubit_no]])


    if wrap ==True:
        circ = circ.to_gate()
        circ.label = "Ry Theta" 

    return circ


def Inspiral_Fixed_Rots(n):
    """Create a circuit for the inspiral f^-7/6
    Args: n: number of qubits
    Returns: Circuit """
    """My theta values?:
    [0.26651054,0.37298941,0.61763383,0.48799182,0.63950672,0.6883351,0.71262227] 
    
    """

    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)
    basic = General_State_Prep([0.26651054,0.37298941,0.61763383,0.48799182,0.63950672,0.6883351,0.71262227])
    circ.append(basic,qr[0:3])
    #This is for m=3
    controls =["000","100","010","110","01","11"]
    theta_m3 =[0.59220374,0.67082012,0.70137334,0.720939,0.73261907,0.74665287]
    m_3=ThetaRotation(circ,qr,controls,3,theta_m3,True)
    circ.append(m_3,[*qr])
    #for the rest of the levels m>3
    for i in range((n-4)):
        controls =["0000","1000","100","10","01","11"]
        thetas=[[0.67107712,0.70335526,0.72139602,0.74093268,0.75778171,0.76536915],[0.72229088,0.7413661,0.75127389,0.76229837,0.7712584,0.77520995],[0.75208093,0.76253209,0.76820956,0.77361685,0.77824197,0.78025934],[0.76825558,0.77373839,0.77667612,0.77944766,0.78179804,0.7828174],[0.77669981,0.77950966,0.78100439,0.7824077,0.78359254,0.78410493]]
        fixed=ThetaRotation(circ,qr,controls,(4+i),thetas[i],True)
        circ.append(fixed,[*qr])
    circ.decompose().draw("mpl",fold=-1)
    circ.save_statevector()
    circ.measure(qr,cla_reg)
    shots = 1000
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    state_vector = result.get_statevector(tqc)
    print("statevector:",state_vector)
    plt.show()


def Waveform_Fixed_Rots(n):
    """[0.69120865] 
    [0.75495778 0.49805573] 
    [0.73051654 0.75904243 0.70661518 0.66738346] 
    [0.72253677 0.74907001 0.75904243 0.70661518 0.66738346] 
    [0.73570614 0.74916687 0.74907001 0.75904243 0.70661518 0.66738346]  """

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
    basic = General_State_Prep([0.74928181,0.73045836,0.72835018,0.71932391,0.73274099,0.74580992,0.75401846])
    
    circ.append(basic,qr[0:3])
    circ.save_statevector(label='v1')
    #This is for m=3
    #circ.save_statevector()
    controls=["000"]
    theta_m3=[0.72]
    m_3=ThetaRotation(circ,qr,controls,3,theta_m3,True)
    circ.append(m_3,[*qr])
    controls =["000","100","10","01","11"]
    #theta_m3=[0.72251893,0.74890938,0.75837712,0.70859427,0.66545498]
    theta_m3 = [0.72434959,0.74059692,0.75007482,0.76163211,0.76743471] 
    m_3=ThetaRotation(circ,qr,controls,3,theta_m3,True)
    circ.append(m_3,[*qr])
    #for the rest of the levels m>3
    circ.save_statevector(label='v2')
    
    for i in range((n-4)):
        controls =["0000","1000","100","10","01","11"]
        #theta_m4=[0.73581314,0.74934639,0.74890938,0.75837712,0.70859427,0.66545498]
        theta_m4=[[0.74066914,0.75075189,0.75726406,0.76502399,0.77230382,0.77575838],[0.75760343,0.76520445,0.76942786,0.77434713,0.77851473,0.78040074],[0.76978013,0.77445391,0.77707962,0.77963685,0.78186772,0.78285331],[0.77710115,0.77969497,0.78110797,0.78245585,0.78361015,0.78411397],[0.78111943,0.78248617,0.78321923,0.77277432,0.78449838,0.7847531]]
        fixed=ThetaRotation(circ,qr,controls,(4+i),theta_m4[i],True)
        circ.append(fixed,[*qr])
    circ.save_statevector(label='v3')
    #circ.save_statevector()    
    '''
    #circ.save_statevector()
    circ.measure(qr,cla_reg)
    backend = qt.Aer.get_backend('statevector_simulator')
    result = qt.execute(circ, backend, shots=500000)
    # Get the statevector from result().
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=50000)
    result = job.result()
    counts = result.get_counts(tqc)
    print("counts:",counts)
    #statevector = result.get_statevector(circ)
    #print(statevector)
    '''
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


def Normalise(values):
    """
    input a string of values
    return: normalised values (add up to 1)
    """
    new = []
    total = sum(values)
    for i in values:
        new.append(i/total)
    print(new)
    return new

def test(n):
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)
    #basic = General_State_Prep([0.69134909,0.75468131,0.498768,0.7304692,0.75837712, 0.70859427, 0.66545498])
    #basic = General_State_Prep([0.74928181,0.73045836,0.72835018,0.71932391,0.73274099,0.74580992,0.75401846,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    basic = General_State_Prep([0.867270505,0.818421983,0.548723411,0.796702905,0.80948195,0.95247867,0.296978124,0.790645952,0.791867976,0.794145441,0.801406291,0.821086042,0.943514493,0.616001409,0.521245153,0.788023101,0.788061049,0.788206058,0.789029605,0.78921198,0.790481467,0.792071168,0.794709994,0.799029946,0.807553099,0.825022378,0.924313151,0.474386805,0.897700349,0.703589962,0.253749431])
    circ.append(basic,qr)
    circ.save_statevector(label='v1')
    circ.decompose().draw("mpl",fold=-1)
    backend = QasmSimulator()
    job = qt.execute(circ, backend, shots=1000)
    result = job.result()
    
    statevector=result.data(0)['v1']
    circ.decompose().draw("mpl",fold=-1)
    plt.show()
    print(statevector)

def amplitudes_for_theta(f):
    """Samples the distributions at the given frequency
    Then normalises them so them squared and summed =1"""
    result =[]
    normalised=[]
    for i in f:
        powered = pow(i,-7/3)
        result.append(powered)
    for x in result:
        norm = np.sqrt(x/sum(result))
        normalised.append(norm)
    return normalised

def amplitudes_7_over_6(f):
    result =[]
    normalised=[]
    for i in f:
        powered = pow(i,-7/6)
        result.append(powered)
    for x in result:
        norm = np.sqrt(x/sum(result))
        normalised.append(norm)
    return normalised

def theta(frequency):
    '''
    Generates thetas for a 
    Args: frequency: A list of all the frequencys you will be calculating theta for  
    '''
    m = 3
    j=0
    #This is specifically for f^-7/6 (so p(i)is f^-7/3)
    #Will need changing for a more complicatated distribution
    amps= amplitudes_for_theta(frequency)
    thetas=[]
    for i in range(m):
        size = int(0+8/((i+1)*2))
        start_f =0
        print("m:",i)
        j=pow(2,i)
        for x in range(j):
            print("j:",x)
            print("size:",size)
            print("start_f:",(start_f))
            print("start_f+size:",(start_f+size))
            print("end:",start_f+size*2)
            upper = sum(amps[start_f:(start_f+size)])
            lesser =sum(amps[start_f:(start_f+size*2)])
            costheta2 = upper/lesser
            costheta = np.sqrt(costheta2)
            theta = np.arccos(costheta)
            
            print("theta:",theta)
            thetas.append(theta*2)
            start_f= start_f+(size*2)
    return thetas

if __name__ == "__main__":
    '''with open('test.csv', newline='') as f:
        reader = csv.reader(f)
        data = list((reader))
    values=[]
    for i in data[0]:
        values.append(float(i))
    Normalise(values)'''
    #Inspiral_Fixed_Rots(9)
    #Waveform_Fixed_Rots(9)
    #test(5)
    #thetas = theta()
    #print(thetas)
    frequency = [40,78.75,117.5,156.25,195,233.75,272.5,311.25]
    thetas = theta(frequency)
    print(thetas)
    #thetas =[[0.46257475965672856],[0.5253295465140008, 0.6946960790369119],[0.5927731247686209, 0.702648001759369, 0.732631041305686, 0.7466576155044997]]
    qr= qt.QuantumRegister(size=3,name='q')
    circ = qt.QuantumCircuit(qr)
    #basic =General_State_Prep([0.26651054,0.37298941,0.61763383,0.48799182,0.63950672,0.6883351,0.71262227])
    basic =General_State_Prep(thetas)
    #x_gate_pattern=[]
    #xGateAdd(2,x_gate_pattern,circ,qr,thetas,3)
    circ.append(basic,qr[0:3])
    circ.save_statevector(label='v1')
    backend = QasmSimulator()
    backend_options = 'statevector'
    job = qt.execute(circ, backend, shots=10000)
    result = job.result()
    statevector=result.data(0)['v1']
    circ.decompose().draw("mpl",fold=-1)
    plt.show()
    print("StateVec:",statevector)
    print("probs:",np.sqrt(statevector.probabilities()))
    amps = amplitudes_for_theta([40,78.75,117.5,156.25,195,233.75,272.5,311.25])
    amps_7_6 = amplitudes_7_over_6([40,78.75,117.5,156.25,195,233.75,272.5,311.25])
    print(amps)
    plt.plot(frequency,np.sqrt(statevector.probabilities()),color='k',label='Statevector')
    plt.plot(frequency,amps,color='r',label='Amps')
    plt.plot(frequency,amps_7_6,color='b',label='Amps')
    plt.legend()
    plt.show()
    test=[1,2,3,4]
    print(sum(test[0:4]))