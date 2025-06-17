import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
import matplotlib.pyplot as plt
import numpy as np


#This is Roselyns code
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
    """
    Generates a circuit made of Rotational Y gates which have controls based on predetermined
    conditions.
    :param circ: quantum circuit
    :param qr: quantum register
    :param condition: list of str what conditions should it have
    :param qubit_no: which qubit
    :param theta:list of corresponding theta values
    :param wrap: boolean wraps gates together into a single gate
    :return: quantum circuit
    """

    size_qr = len(qr)

    #Wraps up the circuit into a single gate
    if wrap == True:
        qr =qt.QuantumRegister(size_qr,"qr")
        circ = qt.QuantumCircuit(qr)
    
    for count,i in enumerate(condition):
        start = len(i)
        rotation_gate = RYGate(theta[count]).control(len(i),ctrl_state=i)
        #[*qr[num_qubits-i-1:][::-1]
        #TODO The 8 here is currently hard coded for n=9 qubits, this could be easily changed by passing through n as a parameter
        circ.append(rotation_gate,[*qr[(8-len(i)+1):],qr[qubit_no]])
        #circ.append(rotation_gate,[*qr[9-3:][::-1]])

    if wrap ==True:
        circ = circ.to_gate()
        circ.label = "Ry Theta" 

    return circ

def amplitudes_for_theta(f):
    """
    Samples the function f^{-7/3} at the given frequency
    Then normalises them so them squared and summed =1
    :param f: list of frequency values
    :return: list of normalised amplitudes
    """
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
    """
    Generates thetas for a list of frequencies
    :param frequency:
    :param m: A list of all the frequencies you will be calculating theta for
    :return: the number of qubits in the circ
    This is specifically for f^-7/6 (so p(i)is f^-7/3)
    TODO Generalise the function for any function
    """
    j=0

    #Gets a list of normalised amplitudes for all the frequencies
    amps= amplitudes_for_theta(frequency)
    thetas=[]
    gsp_theta=[]
    for i in range(m):
        j=pow(2,i)
        #For GSP we calculate all possible theta values
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
        #For our reduced circuit we only calculate theta values at specific starting points (indexs)
        elif i>=3:
            #Our intermediate level
            if i==3:
                thetas.append(gsp_theta)
                '''Starting points frequency: 40,78.75,117,156.25,195,272.5
                These indexs are hard coded for when n=9 would need to be changed for different points
                Also these directly correspond with the set up of the circuit 
                so would also need to be changed if the circuit changes
                TODO Generalise this so we send through index'''
                index=[0,64,128,192,256,384]
                #index=[0,4,8,12,16,24]
            #Final circuit structure
            else:
                '''Set starting points for all the
                40, 59.375, 78.75, 117.5, 195, 272.5
                0,32,64,128,256,384'''
                #TODO Gereralise this function
                index=[0,32,64,128,192,256,384]
                #index=[0,2,4,8,16,24]
            temp_theta =[]
            #Calculates theta for a given index
            for current_index in index:
                temp_start=current_index
                place = current_index/pow(2,m-i)
                increment = int(pow(2,m)/pow(2,i))
                mid =int((place+0.5)*pow(2,m-i))
                end=int((place+1)*pow(2,m-i))
                upper = sum(amps[current_index:mid])
                lesser =sum(amps[current_index:end])
                costheta2 = upper/lesser
                costheta = np.sqrt(costheta2)
                theta = np.arccos(costheta)
                temp_theta.append(theta*2)
            thetas.append(temp_theta)
    if i<3:
        thetas.append(gsp_theta)
    return thetas

def amplitudes_7_over_6(f):
    """
    Calcultes the amplitudes for the function f^{-7/6}
    then normalises them
    :param f: Frequency range
    :return: list of normalised amplitudes
    """
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
    """
    Caluclates the fidelity of two systems.
    :param expected_amps: list of the amplitudes from the discritsed function
    :param measured_amps: list of amplitudes taken from the statevector of the quantum circuit
    :return: the fidelity (float)
    """
    current =0
    for count,i in enumerate(expected_amps):
        current = current+(i*measured_amps[count])
    
    fidelity = current*current
    return fidelity

def Inspiral_Fixed_Rots(n,PLOT=True):
    """
    Create a quantum circuit for the inspiral f^-7/6
    :param n: number of qubits
    :param PLOT: Boolean. Should you display the outcome?
    :return: None
    """

    #Sets up the circuit with a  quantum register n large
    #And a classic Register for measuring
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)

    #Generate theta values for the GSP section
    #A list of frequency values evenly spaced, the amount is based on n
    frequency = np.linspace(40,350,num=pow(2,n),endpoint=False)
    #Generates a two-dimensional array, 0 position will contain all the thetas for GSP
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
        controls =["0000","0001","001","010","011","10","11"]
        fixed=ThetaRotation(circ,qr,controls,8-(4+i),thetas[i+2],True)
        circ.append(fixed,[*qr])
    #circ.decompose().draw("mpl",fold=-1)

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
    #plt.show()

    #Generates the expected amplitudes
    amps = amplitudes_for_theta(frequency)
    amps_7_6 = amplitudes_7_over_6(frequency)

    #Plots the results in comparision to the actual amplitudes
    if PLOT:
        rc_results = np.sqrt(state_vector.probabilities())
        result_fig = plt.figure()
        result_fig.set_figwidth(8)
        result_fig.set_figheight(6)
        plt.plot(frequency,np.sqrt(state_vector.probabilities()),color='k',label='Statevector',ms = 2)
        plt.plot(frequency,amps_7_6,color='r',label='Amps -7/6',linestyle='dashed')
        plt.legend()
        plt.xlabel('f (Hz)')
        plt.ylabel('A')
        result_fig.savefig("../Images/Inspiral_Amplitudes.png",dpi=1000,bbox_inches='tight')
        plt.show()

        #Count the number of gates for the circuit
        print(dict(circ.decompose().count_ops()))
        print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))

        #Plot the box plot to see the errors
        difference = rc_results - amps_7_6
        fig = plt.figure()
        fig.set_figwidth(8)
        fig.set_figheight(2)
        plt.boxplot(abs(difference),vert=False,widths=1)
        plt.ylim(0,2)
        plt.yticks([])
        plt.xlabel("Error")
        #plt.title("Boxplot showing the errors between the datapoints encoded by the reduced circuit\n and the intended amplitudes.")
        fig.savefig('../Images/Inspiral_boxplot.png', dpi=1000,bbox_inches='tight')
        plt.show()

    print("Fidelity: ",Fidelity(amps_7_6,np.sqrt(state_vector.probabilities())))

def theata_Fixed(frequency,m):
    """
    Generates the corresponding theta values for the fixed circuit
    :param frequency:
    :param m: the levels we are going to
    :return: list of all the theta values for the circuit
    """
    #TODO Combine with Theta function
    j=0
    amps= amplitudes_for_theta(frequency)
    thetas=[]
    gsp_theta=[]
    for i in range(m):
        j=pow(2,i)
        #Calculates the Thetas for GSP
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
            #Calculates a single theta for all remaining levels
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

def plot_bounds(x_data,y_data,n):
    """
    Displays G-R Plots for our distribution
    args:
    x_data: our frequency data
    y_data
    n: the number of qubits

    return:
    None

    """

    fmin=0
    fmax=512
    for m in range(n):
        plt.plot(x_data,y_data,color='blue',label='Inspiral')
        print("level:",m)
        temp_fmin = fmin
        freq_increment = (fmax-fmin)/2**(m+1)
        for j in range(2**m):
            print("j:",j)
            print("temp_fmin:",temp_fmin)
            freq_end = temp_fmin +freq_increment
            plt.axvline(x = 195, color = '#635dff')
            plt.fill_between(x_data[int(temp_fmin):int(freq_end)], 0,y_data[int(temp_fmin):int(freq_end)], color='#dcdcfd')
            temp_fmin = freq_end+freq_increment
        plt.legend(["This is my legend"], fontsize="x-large")
        plt.show()


if __name__ == "__main__":
    #Correct one
    Inspiral_Fixed_Rots(9)
    #Generate thetas
    

    #Fixed_Rot_Old(9)
    #test(6)
    '''frequency = np.linspace(40,350,num=pow(2,9),endpoint=False)
    amps_7_6 = amplitudes_7_over_6(frequency)
    plot_bounds(frequency,amps_7_6,2)'''

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
    
