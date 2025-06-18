import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
import matplotlib.pyplot as plt
import numpy as np
import circ_const as CC


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

def theta(frequency,m,circ_type):
    """
    Generates thetas for a list of frequencies
    :param frequency:
    :param m: A list of all the frequencies you will be calculating theta for
    :return: the number of qubits in the circ
    This is specifically for f^-7/6 (so p(i)is f^-7/3)
    """
    #Gets a list of normalised amplitudes for all the frequencies
    '''Starting points frequency: 40,78.75,117,156.25,195,272.5
    These indexs are hard coded for when n=9 would need to be changed for different points'''
    amps= amplitudes_for_theta(frequency)
    if circ_type == "v0":
        thetas= [CC.gsp_thetas(amps, 3, 9)]
        index=[0,64,128,192,256,384]
        thetas.append(CC.indexed_theta(amps,3,9,index))
        index=[0,32,64,128,192,256,384]
        for i in range(5):
            thetas.append(CC.indexed_theta(amps,4+i,9,index))
    elif circ_type =="v1":
        thetas= [CC.gsp_thetas(amps, 3, 9)]
        index=[0,64,128,192,256,320,384]
        thetas.append(CC.indexed_theta(amps,3,9,index))
        index=[0,32,64,96,128,160,192,256,320,384]
        thetas.append(CC.indexed_theta(amps,4,9,index))
        index=[0,16,32,48,64,96,128,160,192,256,320,384]
        for i in range(4):
            thetas.append(CC.indexed_theta(amps,5+i,9,index))

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


def Inspiral_Fixed_Rots(n,gsp,circ_type,PLOT=True):
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
    #The following positions will be containing the thetas for each level up to n-1
    thetas = theta(frequency,n,circ_type)

    #Does the basic GSP algorithm for the first 3 levels
    #m=0 ->m=2
    basic = CC.General_State_Prep(thetas[0])
    circ.append(basic,qr[n-gsp:n])
    if circ_type == 'v0':
        #Currently we have a fixed circuit, the controls are defined here:
        controls =["000","001","010","011","10","11"]
        m_3=CC.ThetaRotation(circ,qr,controls,5,thetas[1],n,True)
        circ.append(m_3,[*qr])
        #for the rest of the levels m>3
        for i in range((n-4)):
            controls =["0000","0001","001","010","011","10","11"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(4+i),thetas[i+2],n,True)
            circ.append(fixed,[*qr])
    elif circ_type == 'v1':
        controls =["000","001","010","011","100","101","11"]
        m_3=CC.ThetaRotation(circ,qr,controls,5,thetas[1],n,True)
        circ.append(m_3,[*qr])
        #for the rest of the levels m>3
        controls =["0000","0001","0010","0011","0100","0101","011","100","101","11"]
        fixed=CC.ThetaRotation(circ,qr,controls,4,thetas[2],n,True)
        circ.append(fixed,[*qr])
        for i in range((n-5)):
            controls =["00000","00001","00010","00011","0010","0011","0100","0101","011","100","101","11"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(5+i),thetas[i+3],n,True)
            circ.append(fixed,[*qr])


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

    print("Fidelity: ",CC.Fidelity(amps_7_6,np.sqrt(state_vector.probabilities())))

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
    Inspiral_Fixed_Rots(9,3,"v0")
    Inspiral_Fixed_Rots(9,3,"v1")
    #Generate thetas
