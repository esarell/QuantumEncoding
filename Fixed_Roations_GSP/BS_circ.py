import qiskit as qt

from numpy.core.defchararray import lower, index
from qiskit.circuit.library.standard_gates import RYGate
from qiskit.providers.aer import QasmSimulator
import matplotlib.pyplot as plt
import numpy as np
import circ_const as CC

def Black_Scholes(x,K,s):
    #Genrates the distribution
    y = []
    test = np.log(K * s)
    print(-test)
    for i in x:
        if -test <= i < 0:
            y.append(K - (np.exp(-i) / s))
        elif 0 < i <= test:
            y.append(K - (np.exp(i) / s))
        else:
            print("inf")
            y.append(0.0000001)
    return y

def BS_angle_Gen(x,y,m,circ_type):
    j=0
    #Gets a list of normalised amplitudes for all the frequencies
    amps=CC.prob_normalised(y)
    thetas=[]
    if circ_type == "v1":
        thetas=[CC.gsp_thetas(amps, 3, 5)]
        index_m = [0,4,8,16,24,28]
        thetas.append(CC.indexed_theta(amps,3,5,index_m))
        index_m = [0,2,4,8,16,24,28,30]
        thetas.append(CC.indexed_theta(amps,4,5,index_m))
    elif circ_type =="v2":
        thetas=[CC.gsp_thetas(amps, 3, 5)]
        index_m = [0,4,8,16,20,24,28]
        thetas.append(CC.indexed_theta(amps,3,5,index_m))
        index_m = [0,2,4,6,8,16,20,24,26,28,30]
        thetas.append(CC.indexed_theta(amps,4,5,index_m))
    elif circ_type ==   "test":
        thetas=[CC.gsp_thetas(amps, 3, 5)]
        index_m = [0,4,8,12,16,20,24,28]
        thetas.append(CC.indexed_theta(amps,3,5,index_m))
        index_m = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
        thetas.append(CC.indexed_theta(amps,4,5,index_m))
    elif circ_type == "v3":
        thetas=[CC.gsp_thetas(amps, 3, 9)]
        index_m3 =[0,64,128,256,384,448]
        thetas.append(CC.indexed_theta(amps,3,9,index_m3))
        index_m4 =[0,32,64,128,256,384,448,480]
        thetas.append(CC.indexed_theta(amps,4,9,index_m4))
        index_m5 =[0,16,32,64,128,256,384,448,480,496]
        for i in range(4):
            thetas.append(CC.indexed_theta(amps,5+i,9,index_m5))
    elif circ_type =="v4":
        thetas=[CC.gsp_thetas(amps, 4, 9)]
        index_m4 =[0,32,64,96,128,192,256,320,384,416,448,480]
        thetas.append(CC.indexed_theta(amps,4,9,index_m4))
        index_m5 =[0,16,32,48,64,96,128,192,256,320,384,416,448,464,480,496]
        thetas.append(CC.indexed_theta(amps,5,9,index_m5))
        index_m6 =[0,8,16,24,32,48,64,96,128,192,256,320,384,416,448,464,480,488,496,504]
        thetas.append(CC.indexed_theta(amps,6,9,index_m6))
        index_m7 =[0,4,8,16,24,32,48,64,96,128,192,256,320,384,416,448,464,480,488,496,504,508]
        for i in range(2):
            thetas.append(CC.indexed_theta(amps,7+i,9,index_m7))
    elif circ_type =="v5":
        thetas=[CC.gsp_thetas(amps, 3, 9)]
        index=[0,64,128,256,384,448]
        thetas.append(CC.indexed_theta(amps,3,9,index))
        index=[0,32,64,128,256,384,448,480]
        thetas.append(CC.indexed_theta(amps,4,9,index))
        index=[0,16,32,64,128,256,384,448,480,496]
        for i in range(4):
            thetas.append(CC.indexed_theta(amps,5+i,9,index))
    elif circ_type =="v6":
        thetas=[CC.gsp_thetas(amps, 3, 9)]
        index=[0,64,128,256,384,448]
        thetas.append(CC.indexed_theta(amps,3,9,index))
        index=[0,32,64,128,256,384,448,480]
        for i in range(5):
            thetas.append(CC.indexed_theta(amps,4+i,9,index))
    elif circ_type =="v7":
        thetas=[CC.gsp_thetas(amps, 3, 9)]
        index=[0,64,128,256,384,448]
        for i in range(6):
            thetas.append(CC.indexed_theta(amps,3+i,9,index))
    elif circ_type =="v8":
        thetas=[CC.gsp_thetas(amps, 3, 12)]
        index=[0,512,1024,2048,3072,3584]
        thetas.append(CC.indexed_theta(amps,3,12,index))
        index=[0,256,512,1024,2048,3072,3584,3840]
        thetas.append(CC.indexed_theta(amps,4,12,index))
        index=[0,128,256,512,1024,2048,3072,3584,3840,3968]
        for i in range(7):
            thetas.append(CC.indexed_theta(amps,5+i,12,index))

    return thetas

def BS_circ(n,circ_type,Plot=False):
    """
    Create a quantum circuit for the Black-Schols Problem
    :param n: number of qubits
    :return:
    """
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)

    #Generate theta values for the GSP section
    #A list of frequency values evenly spaced, the amount is based on n
    x_values = np.linspace(-8,8,num=pow(2,n),endpoint=True)
    print("xlen:",len(x_values))
    y_values = Black_Scholes(x_values,45,(45*3))
    theta_vals = BS_angle_Gen(x_values,y_values,n,circ_type)

    if circ_type == "er":
        gsp = CC.General_State_Prep(theta_vals[0])
        circ.append(gsp,qr[2:n])
        m_3_controls =["000","001","01","10","110","111"]
        m_3=CC.ThetaRotation(circ,qr,m_3_controls,1,theta_vals[1],5,True)
        circ.append(m_3,[*qr])
        m_4_controls =["0000","0001","001","01","10","110","1110","1111"]
        m_4=CC.ThetaRotation(circ,qr,m_4_controls,0,theta_vals[2],5,True)
        circ.append(m_4,[*qr])
        '''
        for i in range((n-5)):
            controls =["00000","00001","0001","001","01","10","110","1110","11110","11111"]
            fixed=CC.ThetaRotation(circ,qr,controls,7-(5+i),theta_vals[i+3],True)
            circ.append(fixed,[*qr])'''
    elif circ_type =="v1":
        gsp = CC.General_State_Prep(theta_vals[0])
        circ.append(gsp,qr[n-3:n])
        controls = ["000","001","01","10","110","111"]
        m_3=CC.ThetaRotation(circ,qr,controls,1,theta_vals[1],5,True)
        circ.append(m_3,[*qr])
        controls = ["0000","0001","001","01","10","110","1110","1111"]
        m_4=CC.ThetaRotation(circ,qr,controls,0,theta_vals[2],5,True)
        circ.append(m_4,[*qr])
    elif circ_type =="v2":
        print("here")
        gsp = CC.General_State_Prep(theta_vals[0])
        circ.append(gsp,qr[n-3:n])
        controls = ["000","001","01","100","101","110","111"]
        m_3=CC.ThetaRotation(circ,qr,controls,1,theta_vals[1],5,True)
        circ.append(m_3,[*qr])
        controls = ["0000","0001","0010","0011","01","100","101","1100","1101","1110","1111"]
        m_4=CC.ThetaRotation(circ,qr,controls,0,theta_vals[2],5,True)
        circ.append(m_4,[*qr])
    elif circ_type == "v3":
        gsp = CC.General_State_Prep(theta_vals[0])
        circ.append(gsp,qr[n-3:n])
        controls = ["000","001","01","10","110","111"]
        m_3=CC.ThetaRotation(circ,qr,controls,5,theta_vals[1],9,True)
        circ.append(m_3,[*qr])
        controls = ["0000","0001","001","01","10","110","1110","1111"]
        m_4=CC.ThetaRotation(circ,qr,controls,4,theta_vals[2],9,True)
        circ.append(m_4,[*qr])
        for i in range(n-5):
            controls = ["00000","00001","0001","001","01","10","110","1110","11110","11111"]
            fixed_circ = CC.ThetaRotation(circ,qr,controls,8-(5+i),theta_vals[i+3],n,True)
            circ.append(fixed_circ,[*qr])
    elif circ_type == "v4":
        gsp = CC.General_State_Prep(theta_vals[0])
        circ.append(gsp,qr[n-4:n])
        controls = ["0000","0001","0010","0011","010","011","100","101","1100","1101","1110","1111"]
        m_4=CC.ThetaRotation(circ,qr,controls,4,theta_vals[1],n,True)
        circ.append(m_4,[*qr])
        controls = ["00000","00001","00010","00011","0010","0011","010","011","100","101","1100","1101","11100","11101","11110","11111"]
        m_5 = CC.ThetaRotation(circ,qr,controls,3,theta_vals[2],n,True)
        circ.append(m_5,[*qr])
        controls = ["000000","000001","000010","000011","00010","00011","0010","0011","010","011","100","101","1100","1101","11100","11101","111100","111101","111110","111111"]
        m_6 = CC.ThetaRotation(circ,qr,controls,2,theta_vals[3],n,True)
        circ.append(m_6,[*qr])
        controls = ["0000000","0000001","000001","000010","000011","00010","00011","0010","0011","010","011","100","101","1100","1101","11100","11101","111100","111101","111110","1111110","1111111"]
        m_7 = CC.ThetaRotation(circ,qr,controls,1,theta_vals[4],n,True)
        circ.append(m_7,[*qr])
        m_8 = CC.ThetaRotation(circ,qr,controls,0,theta_vals[5],n,True)
        circ.append(m_8,[*qr])
    elif circ_type =="v5":
        basic = CC.General_State_Prep(theta_vals[0])
        circ.append(basic,qr[n-3:n])
        controls =["000","001","01","10","110","111"]
        m_3=CC.ThetaRotation(circ,qr,controls,5,theta_vals[1],n,True)
        circ.append(m_3,[*qr])
        controls =["0000","0001","001","01","10","110","1110","1111"]
        m_3=CC.ThetaRotation(circ,qr,controls,4,theta_vals[2],n,True)
        circ.append(m_3,[*qr])
        for i in range((n-5)):
            controls =["00000","00001","0001","001","01","10","110","1110","11110","11111"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(5+i),theta_vals[i+3],n,True)
            circ.append(fixed,[*qr])
    elif circ_type =="v6":
        basic = CC.General_State_Prep(theta_vals[0])
        circ.append(basic,qr[n-3:n])
        controls =["000","001","01","10","110","111"]
        m_3=CC.ThetaRotation(circ,qr,controls,5,theta_vals[1],n,True)
        circ.append(m_3,[*qr])
        controls =["0000","0001","001","01","10","110","1110","1111"]
        for i in range((n-4)):
            fixed=CC.ThetaRotation(circ,qr,controls,8-(4+i),theta_vals[i+2],n,True)
            circ.append(fixed,[*qr])
    elif circ_type == "v7":
        basic = CC.General_State_Prep(theta_vals[0])
        circ.append(basic,qr[n-3:n])
        controls =["000","001","01","10","110","111"]
        for i in range((n-3)):
            fixed=CC.ThetaRotation(circ,qr,controls,8-(3+i),theta_vals[i+1],n,True)
            circ.append(fixed,[*qr])
    elif circ_type =="v8":
        basic = CC.General_State_Prep(theta_vals[0])
        circ.append(basic,qr[n-3:n])
        controls =["000","001","01","10","110","111"]
        m_3=CC.ThetaRotation(circ,qr,controls,8,theta_vals[1],n,True)
        circ.append(m_3,[*qr])
        controls =["0000","0001","001","01","10","110","1110","1111"]
        m_3=CC.ThetaRotation(circ,qr,controls,7,theta_vals[2],n,True)
        circ.append(m_3,[*qr])
        for i in range((n-5)):
            controls =["00000","00001","0001","001","01","10","110","1110","11110","11111"]
            fixed=CC.ThetaRotation(circ,qr,controls,11-(5+i),theta_vals[i+3],n,True)
            circ.append(fixed,[*qr])

    circ.save_statevector(label='end')
    backend = QasmSimulator()
    backend_options = 'statevector'
    job = qt.execute(circ, backend, shots=1000)
    result = job.result()
    statevector3=result.data(0)['end']
    circ.decompose().draw("mpl",fold=-1)
    plt.show()

    #amps = waveform_amps()
    temp_amps= y_values
    amps=CC.amp_normalised(temp_amps)
    #amps_7_6 = amplitudes_7_over_6(frequency)

    if Plot:
        rc_results = np.sqrt(statevector3.probabilities())
        result_fig = plt.figure()
        result_fig.set_figwidth(8)
        result_fig.set_figheight(6)
        plt.plot(x_values,rc_results,color='k',label='Statevector',ms = 2)
        plt.plot(x_values,amps,color='r',label='Black-Scholes',linestyle='dashed')

        '''plt.axvline(x = -8, color = 'b',label='Controls')
        plt.axvline(x = -7.5, color = 'b')
        plt.axvline(x = -7, color = 'b',)
        plt.axvline(x = -6, color = 'b',)
        plt.axvline(x = -4, color = 'b',)
        plt.axvline(x = 0, color = 'b',)
        plt.axvline(x = 4, color = 'b')
        plt.axvline(x = 6, color = 'b')
        plt.axvline(x = 7, color = 'b')
        plt.axvline(x = 7.5, color = 'b')'''

        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        result_fig.savefig("../Images/Black_Scholes_Amplitudes_"+circ_type+".png",dpi=1000,bbox_inches='tight')
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
        fig.savefig('../Images/Black_Scholes_boxplot_'+circ_type+'.png', dpi=1000,bbox_inches='tight')
        plt.show()


    print("Fidelity: ",CC.Fidelity(amps,np.sqrt(statevector3.probabilities())))
    print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))

def BS_GR_Circ(n):
    """
    Black-Scholes encoding using the Grover-Rudolph Algorithm
    :param n: number of qubits
    :return:
    """
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)
    #Generate theta values for the GSP section
    #A list of frequency values evenly spaced, the amount is based on n
    x_values = np.linspace(-8,8,num=pow(2,n),endpoint=True)
    print("xlen:",len(x_values))
    y_values = Black_Scholes(x_values,45,(45*3))
    amps=CC.prob_normalised(y_values)
    #Change this if you want more than 9 qubits
    theta_vals =[CC.gsp_thetas(amps,9 , 9)]
    amps=CC.amp_normalised(y_values)
    gsp = CC.General_State_Prep(theta_vals[0])
    print("here1")
    circ.append(gsp,qr[:])
    circ.save_statevector(label='end')
    backend = QasmSimulator()
    backend_options = 'statevector'
    job = qt.execute(circ, backend, shots=1000)
    result = job.result()
    statevector3=result.data(0)['end']
    #circ.decompose().draw("mpl",fold=-1)
    #plt.show()
    rc_results = np.sqrt(statevector3.probabilities())
    result_fig = plt.figure()
    result_fig.set_figwidth(8)
    result_fig.set_figheight(6)
    plt.plot(x_values,np.sqrt(statevector3.probabilities()),color='k',label='Statevector',ms = 2)
    plt.plot(x_values,amps,color='r',label='Goal',linestyle='dashed')
    plt.legend()
    plt.xlabel('X')
    plt.ylabel('Amplitude')
    plt.show()


    print("Fidelity: ",CC.Fidelity(amps,np.sqrt(statevector3.probabilities())))
    print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))



if __name__ == "__main__":
    #BS_circ(5,"v1",Plot=True)
    #BS_circ(5,"v2",Plot=True)
    #BS_circ(12,"v8",Plot=True)
    BS_GR_Circ(9)

