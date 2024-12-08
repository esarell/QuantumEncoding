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

    #Wraps up the circuit into a single gate
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

def theta(frequency,m):
    '''
    Generates thetas for a list of frequencys
    Args: frequency: A list of all the frequencys you will be calculating theta for
    m: the number of qubits in the circ
    '''
    j=0
    #This is specifically for f^-7/6 (so p(i)is f^-7/3)
    #Will need changing for a more complicatated distribution
    #Gets a list of normalised amplitudes for all the frequencies 
    amps= amplitudes_for_theta(frequency)
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
                '''Starting points frequency: 40,78.75,117,156.25,195,272.5
                #These indexs are hard coded for when n=9 would need to be changed for different points
                #Also these directly correspond with the set up of the circuit so would also need to be changed if the circuit changes'''
                index=[0,64,128,192,256,384]
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



def Inspiral_Fixed_Rots(n):
    """Create a circuit for the inspiral f^-7/6
    Args: n: number of qubits """

    #Creates the circuit, quantum register n large
    #And a classic Register for measuring
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)

    #Generate theta values for the GSP section
    #A list of frequency values evenly spaced, the amount is based on n
    frequency = np.linspace(40,350,num=pow(2,n),endpoint=False)
    #Generates a two dimensional array, 0 position will contain all the thetas for GSP
    #The following positions will containing the theatas for each level up to n-1
    thetas = theta(frequency,n)

    #Does the basic GSP algorithm for the first 3 levels
    #m=0 ->m=2
    basic = General_State_Prep(thetas[0])
    circ.append(basic,qr[n-3:n])
    #This is for m=3
    #Currently we have a fixed circuit, the controls are defined here:
    controls =["000","001","010","011","10","11"]
    m_3=ThetaRotation(circ,qr,controls,5,thetas[1],True)
    circ.append(m_3,[*qr])
    #for the rest of the levels m>3
    for i in range((n-4)):
        controls =["0000","0001","001","01","10","11"]
        fixed=ThetaRotation(circ,qr,controls,8-(4+i),thetas[i+2],True)
        circ.append(fixed,[*qr])
    circ.decompose().draw("mpl",fold=-1)

    #Get the statevector at the current point in the circuit
    circ.save_statevector()
    #Perform a measurement
    circ.measure(qr,cla_reg)
    shots = 1000
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    state_vector = result.get_statevector(tqc)
    #print("statevector:",state_vector)
    plt.show()

    amps = amplitudes_for_theta(frequency)
    amps_7_6 = amplitudes_7_over_6(frequency)
    plt.plot(frequency,np.sqrt(state_vector.probabilities()),color='k',label='Statevector',marker = 'o',ms = 2,linestyle='None')
    #plt.plot(frequency,amps,color='r',label='Amps -7/3')
    plt.plot(frequency,amps_7_6,color='b',label='Amps -7/6')
    plt.legend()
    plt.xlabel('f (Hz)')
    plt.ylabel('A')
    plt.show()
    print("Fidelity: ",Fidelity(amps_7_6,np.sqrt(state_vector.probabilities())))
    print(dict(circ.decompose().count_ops()))
    print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))

def theata_Fixed(frequency,m):
    j=0
    amps= amplitudes_for_theta(frequency)
    thetas=[]
    gsp_theta=[]
    for i in range(m):
        j=pow(2,i)
        if i<6:
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
        else:
            if i ==6:
                thetas.append(gsp_theta)
            current_index=0
            place = current_index/pow(2,m-i)
            increment = int(pow(2,m)/pow(2,i))
            mid =int((place+0.5)*pow(2,m-i))
            end=int((place+1)*pow(2,m-i))
            upper = sum(amps[current_index:mid])
            lesser =sum(amps[current_index:end])
            costheta2 = upper/lesser
            costheta = np.sqrt(costheta2)
            theta = np.arccos(costheta) 
            thetas.append(theta*2)
    return thetas


def Fixed_Rot_Old(n):
    ''' Fixed Rotations go up to level k then go a single rotation
    '''
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)
    #had to do a reduced frequency range for fixed rotations
    frequency = np.linspace(40,170,num=pow(2,n),endpoint=False)
    #k=is 6
    thetas = theata_Fixed(frequency,n)
    print(thetas)
    basic = General_State_Prep(thetas[0])
    circ.append(basic,qr[n-6:n])
    print(thetas[2])
    for i in range((n-6)):
        rotation_gate = RYGate(thetas[i+1])
        print(n-7-i)
        circ.append(rotation_gate,[qr[n-7-i]])
    circ.save_statevector()

    circ.measure(qr,cla_reg)
    shots = 1000
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    state_vector = result.get_statevector(tqc)
    #print("statevector:",state_vector)
    circ.decompose().draw("mpl",fold=-1)
    plt.show()
    amps = amplitudes_for_theta(frequency)
    amps_7_6 = amplitudes_7_over_6(frequency)
    plt.plot(frequency,np.sqrt(state_vector.probabilities()),color='k',label='Fixed')
    #plt.plot(frequency,amps,color='r',label='Amps -7/3')
    plt.plot(frequency,amps_7_6,color='b',label='Amps -7/6')
    plt.legend()
    plt.xlabel('f (Hz)')
    plt.ylabel('A')
    plt.show()
    print("Fidelity: ",Fidelity(amps_7_6,np.sqrt(state_vector.probabilities())))
    print(dict(circ.decompose().count_ops()))
    print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))

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
    elements = np.linspace(0., 3., 63)
    #basic = General_State_Prep([0.867270505,0.818421983,0.548723411,0.796702905,0.80948195,0.95247867,0.296978124,0.790645952,0.791867976,0.794145441,0.801406291,0.821086042,0.943514493,0.616001409,0.521245153,0.788023101,0.788061049,0.788206058,0.789029605,0.78921198,0.790481467,0.792071168,0.794709994,0.799029946,0.807553099,0.825022378,0.924313151,0.474386805,0.897700349,0.703589962,0.253749431])
    basic = General_State_Prep(elements)
    circ.append(basic,qr)
    '''circ.save_statevector(label='v1')
    circ.decompose().draw("mpl",fold=-1)'''
    backend = QasmSimulator()
    job = qt.execute(circ, backend, shots=1000)
    result = job.result()
    
    #statevector=result.data(0)['v1']
    circ.decompose().draw("mpl",fold=-1)
    plt.show()
    print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))
    #print(statevector)

def plot_bounds(x_data,y_data,n):
    """
    args:
    x_data
    y_data
    n: the level
    """
    fmin=0
    fmax=512
    for m in range(n):
        plt.plot(frequency,amps_7_6,color='blue',label='Inspiral')
        print("level:",m)
        temp_fmin = fmin
        freq_increment = (fmax-fmin)/2**(m+1)
        for j in range(2**m):
            print("j:",j)
            print("temp_fmin:",temp_fmin)
            freq_end = temp_fmin +freq_increment
            plt.axvline(x = 195, color = '#635dff')
            plt.fill_between(frequency[int(temp_fmin):int(freq_end)], 0,amps_7_6[int(temp_fmin):int(freq_end)], color='#dcdcfd')
            temp_fmin = freq_end+freq_increment
        plt.legend(["This is my legend"], fontsize="x-large")
        plt.show()

    


if __name__ == "__main__":
    '''with open('test.csv', newline='') as f:
        reader = csv.reader(f)
        data = list((reader))
    values=[]
    for i in data[0]:
        values.append(float(i))
    Normalise(values)'''
    frequency = np.linspace(40,350,num=pow(2,9),endpoint=False)
    amps_7_6 = amplitudes_7_over_6(frequency)
    print(Fidelity(amps_7_6,amps_7_6))
    Inspiral_Fixed_Rots(9)
    #Fixed_Rot_Old(9)
    #test(6)
    frequency = np.linspace(40,350,num=pow(2,9),endpoint=False)
    amps_7_6 = amplitudes_7_over_6(frequency)
    plot_bounds(frequency,amps_7_6,2)
    '''m=3
    plt.axvline(x = 195, color = '#635dff', label = 'Reduced')
    plt.axvline(x = 117.5, color = '#54d5c2', label = 'Complete')
    plt.axvline(x = 272.5, color = '#54d5c2')

    plt.axvline(x = 78.75, color = '#54d5c2',)
    plt.axvline(x = 156.25, color = '#54d5c2')
    #plt.axvline(x = 233.75, color = 'black')
    plt.axvline(x = 272.5, color = '#635dff')



    plt.fill_between(frequency[0:32], 0,amps_7_6[0:32], color='#dcdcfd',label=r'$\cos^2(\theta)$')
    plt.fill_between(frequency[64:96], 0,amps_7_6[64:96], color='#dcdcfd')
    plt.fill_between(frequency[128:160], 0,amps_7_6[128:160], color='#dcdcfd')
    plt.fill_between(frequency[192:224], 0,amps_7_6[192:224], color='#dcdcfd')
    plt.fill_between(frequency[256:288], 0,amps_7_6[256:288], color='#dcdcfd')
    plt.fill_between(frequency[384:416], 0,amps_7_6[384:416], color='#dcdcfd')'''
    #m=4
    '''plt.axvline(x = 195, color = '#ec3e3e', label = 'Fixed Bounds')
    plt.axvline(x = 272.5, color = '#ec3e3e')

    plt.axvline(x = 117.5, color = '#ec3e3e')
    plt.axvline(x = 272.5, color = '#ec3e3e')
    
    plt.axvline(x = 78.75, color = '#ec3e3e')
    plt.axvline(x = 59.375, color = '#ec3e3e')'''
    
    #plt.axvline(x = 233.75, color = 'black')
    



    '''plt.fill_between(frequency[0:16], 0,amps_7_6[0:16], color='#ffd9d9',label=r'$\cos^2(\theta)$')
    plt.fill_between(frequency[32:48], 0,amps_7_6[32:48], color='#ffd9d9')
    plt.fill_between(frequency[64:80], 0,amps_7_6[64:80], color='#ffd9d9')
    plt.fill_between(frequency[128:144], 0,amps_7_6[128:144], color='#ffd9d9')
    
    plt.fill_between(frequency[256:272], 0,amps_7_6[256:272], color='#ffd9d9')
    plt.fill_between(frequency[384:400], 0,amps_7_6[384:400], color='#ffd9d9')'''
    '''plt.fill_between(frequency[0:8], 0,amps_7_6[0:8], color='#ffd9d9',label=r'$\cos^2(\theta)$')
    plt.fill_between(frequency[32:40], 0,amps_7_6[32:40], color='#ffd9d9')
    plt.fill_between(frequency[64:72], 0,amps_7_6[64:72], color='#ffd9d9')
    plt.fill_between(frequency[128:136], 0,amps_7_6[128:136], color='#ffd9d9')
    
    plt.fill_between(frequency[256:264], 0,amps_7_6[256:264], color='#ffd9d9')
    plt.fill_between(frequency[384:392], 0,amps_7_6[384:392], color='#ffd9d9')'''
    plt.title("Inspiral")
    plt.legend(loc='upper right')
    plt.xlabel('f (Hz)')
    plt.ylabel('A')
    plt.show()
    #Waveform_Fixed_Rots(9)
    #test(5)
    #thetas = theta()
    #print(thetas)
    #frequency = [40,78.75,117.5,156.25,195,233.75,272.5,311.25]
    '''frequency = np.linspace(40,350,num=512,endpoint=False)
    print(frequency)
    thetas = theta(frequency)'''
    '''qr= qt.QuantumRegister(size=3,name='q')
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
    amps = amplitudes_for_theta(frequency)
    amps_7_6 = amplitudes_7_over_6(frequency)
    #plt.plot(frequency,np.sqrt(statevector.probabilities()),color='k',label='Statevector')
    plt.plot(frequency,amps,color='r',label='Amps -7/3')
    plt.plot(frequency,amps_7_6,color='b',label='Amps -7/6')
    plt.legend()
    plt.show()'''

    
