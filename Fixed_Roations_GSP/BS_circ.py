import qiskit as qt
from numpy.core.defchararray import lower
from qiskit.circuit.library.standard_gates import RYGate
from qiskit.providers.aer import QasmSimulator
import matplotlib.pyplot as plt
import numpy as np
import circ_const as CC

def Black_Scholes(x,K,s):
    #Genrates the distribution
    y = []
    test = np.log(K*s)
    print(-test)
    for i in x:
        if -test <= i < 0:
            y.append( K-(np.exp(-i)/s))
        elif 0 < i <= test:
            y.append( K-(np.exp(i)/s))
        else:
            y.append(0.0000001)
    return y

def BS_angle_Gen(x,y,m):
    j=0
    #Gets a list of normalised amplitudes for all the frequencies
    amps=(y)
    print(amps)
    thetas=[]
    gsp_theta=[]
    for i in range(m):
        print(i)
        j=pow(2,i)
        if i<3:
            for x in range(j):
                start=int(x*pow(2,m-i))
                mid = int((x+0.5)*pow(2,m-i))
                end = int((x+1)*pow(2,m-i))
                upper = sum(amps[start:mid])
                lesser =sum(amps[start:end])
                costheta2 = upper/lesser
                costheta = np.sqrt(costheta2)
                theta = np.arccos(costheta)
                gsp_theta.append(theta*2)
            print(j)
        elif i>=3:
            if i == 3:
                thetas.append(gsp_theta)
                #index=[0,32,64,128,192,224]
                index=[0,4,8,16,24,28]
            elif i==4:
                #These indexs are hard coded for when n=9 would need to be changed for different points
                #Also these directly correspond with the set up of the circuit so would also need to be changed if the circuit changes
                #index=[0,16,32,64,128,192,224,240]
                index=[0,2,4,8,16,24,28,30]
            else:
                index=[0,8,16,32,64,128,192,224,240,248]
            temp_theta =[]
            for current_index in index:
                print("current_index",current_index)
                temp_start=current_index
                place = current_index/pow(2,m-i)
                increment = int(pow(2,m)/pow(2,i))
                mid =int((place+0.5)*pow(2,m-i))
                end=int((place+1)*pow(2,m-i))
                print("amps_len:",len(amps))
                upper = sum(amps[current_index:mid])
                print("upper,",upper)
                lesser =sum(amps[current_index:end])
                print("lower",lesser)
                costheta2 = upper/lesser
                costheta = np.sqrt(costheta2)
                theta = np.arccos(costheta)
                temp_theta.append(theta*2)
            thetas.append(temp_theta)
    return thetas

def BS_circ(n,Plot=False):
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
    y_values = Black_Scholes(x_values,45,3)

    x_theta = np.linspace(-8,8,num=pow(2,11),endpoint=True)
    y_theta = Black_Scholes(x_theta,45,3)
    theta_vals = BS_angle_Gen(x_values,y_values,n)
    print('len:',(len(theta_vals)))
    gsp = CC.General_State_Prep(theta_vals[0])
    circ.append(gsp,qr[2:n])

    #This is for m=3
    #circ.save_statevector()
    m_3_controls =["000","001","01","10","110","111"]

    m_3=CC.ThetaRotation(circ,qr,m_3_controls,1,theta_vals[1],True)
    circ.append(m_3,[*qr])

    m_4_controls =["0000","0001","001","01","10","110","1110","1111"]
    m_4=CC.ThetaRotation(circ,qr,m_4_controls,0,theta_vals[2],True)
    circ.append(m_4,[*qr])
    '''
    for i in range((n-5)):
        controls =["00000","00001","0001","001","01","10","110","1110","11110","11111"]
        fixed=CC.ThetaRotation(circ,qr,controls,7-(5+i),theta_vals[i+3],True)
        circ.append(fixed,[*qr])'''

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
        plt.plot(x_values,amps,color='r',label='Amps -7/6',linestyle='dashed')
        plt.legend()
        plt.xlabel('f (Hz)')
        plt.ylabel('A')
        result_fig.savefig("../Images/Black_Scholes_Amplitudes.png",dpi=1000,bbox_inches='tight')
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
        fig.savefig('../Images/Black_Scholes_boxplot.png', dpi=1000,bbox_inches='tight')
        plt.show()


    print("Fidelity: ",CC.Fidelity(amps,np.sqrt(statevector3.probabilities())))
    print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))


if __name__ == "__main__":
    BS_circ(5,Plot=True)