import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
from qiskit.providers.aer import QasmSimulator
import matplotlib.pyplot as plt
import numpy as np
import pickle
import circ_const as CC

def waveform_amps():
    """
    :return:
    """
    with open("waveform_amps.pkl", 'rb') as pickleFile:
        test = pickle.load(pickleFile)
    return test

def waveform_theta(frequency,m,circ_type):
    """

    :param frequency:
    :param m:
    :return:
    """
    #Gets a list of normalised amplitudes for all the frequencies 
    temp_amps= waveform_amps()
    amps=CC.prob_normalised(temp_amps)
    if circ_type == "v0":
        thetas=[CC.gsp_thetas(amps, 4, 9)]
        index =[0,32,64,128,192,256,320,352,384,416,448]
        for i in range(5):
            thetas.append(CC.indexed_theta(amps,4+i,9,index))
    elif circ_type == "v1":
        thetas=[CC.gsp_thetas(amps, 4, 9)]
        index =[0,32,64,96,128,160,192,256,320,352,384,416,448,480]
        thetas.append(CC.indexed_theta(amps,4,9,index))
        index=[0,16,32,48,64,80,96,128,160,192,256,320,336,352,368,384,400,416,432,448,480]
        for i in range(4):
            thetas.append(CC.indexed_theta(amps,5+i,9,index))
    return thetas

def Waveform_Fixed_Rots(n,circ_type,Plot=False):
    """
    Create a quantum circuit for the whole gravitational waveform
    Waveform can be seen in the LinearSplineApprox.py file
    :param n: number of qubits
    :return:
    """
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)

    #Generate theta values for the GSP section
    #A list of frequency values evenly spaced, the amount is based on n
    frequency = np.linspace(40,350,num=pow(2,n),endpoint=False)
    angles = waveform_theta(frequency,n,circ_type)

    if circ_type == "v0":
        GSP = CC.General_State_Prep(angles[0])
        circ.append(GSP,qr[n-4:n])
        #for the rest of the levels m>3
        for i in range((n-4)):
            controls =["0000","0001","001","010","011","100","1010","1011","1100","1101","111"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(4+i),angles[i+1],n,True)
            circ.append(fixed,[*qr])
    elif circ_type == "v1":
        GSP = CC.General_State_Prep(angles[0])
        circ.append(GSP,qr[n-4:n])
        controls =["0000","0001","0010","0011","0100","0101","011","100","1010","1011","1100","1101","1110","1111"]
        fixed=CC.ThetaRotation(circ,qr,controls,4,angles[1],n,True)
        circ.append(fixed,[*qr])
        for i in range(n-5):
            controls = ["00000","00001","00010","00011","00100","00101","0011","0100","0101","011","100","10100","10101","10110","10111","11000","11001","11010","11011","1110","1111"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(5+i),angles[i+2],n,True)
            circ.append(fixed,[*qr])

    circ.save_statevector(label='end')
    backend = QasmSimulator()
    backend_options = 'statevector'
    job = qt.execute(circ, backend, shots=1000)
    result = job.result()
    statevector3=result.data(0)['end']
    circ.decompose().draw("mpl",fold=-1)
    plt.show()

    temp_amps= waveform_amps()
    amps=CC.amp_normalised(temp_amps)

    if Plot:
        rc_results = np.sqrt(statevector3.probabilities())
        result_fig = plt.figure()
        result_fig.set_figwidth(8)
        result_fig.set_figheight(6)
        plt.plot(frequency,rc_results,color='k',label='Statevector',ms = 2)
        plt.plot(frequency,amps,color='r',label='Amps -7/6',linestyle='dashed')
        plt.legend()
        plt.xlabel('f (Hz)')
        plt.ylabel('A')
        result_fig.savefig("../Images/Waveform_Amplitudes.png",dpi=1000,bbox_inches='tight')
        plt.show()

        #Plot the box plot to see the errors
        difference = rc_results - amps
        fig = plt.figure()
        fig.set_figwidth(8)
        fig.set_figheight(2)
        plt.boxplot(abs(difference),vert=False,widths=1)
        plt.ylim(0,2)
        plt.yticks([])
        plt.xlabel("Error")
        #plt.title("Boxplot showing the errors between the datapoints encoded by the reduced circuit\n and the intended amplitudes.")
        fig.savefig('../Images/Waveform_boxplot.png', dpi=1000,bbox_inches='tight')
        plt.show()

    print("Fidelity: ",CC.Fidelity(amps,np.sqrt(statevector3.probabilities())))
    print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))

if __name__ == "__main__":
    Waveform_Fixed_Rots(9,"v0",Plot=True)
    Waveform_Fixed_Rots(9,"v1",Plot=True)
