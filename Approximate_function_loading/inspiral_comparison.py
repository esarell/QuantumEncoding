'''This is to generated the reduced circuit for a small frequency range
40-170 Hz'''

import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
import matplotlib.pyplot as plt
import numpy as np
import math
import circ_const as CC
import inspiral_afl as IAFL

def gen_angles(frequency,n,circ_type):
    amps= IAFL.amplitudes_for_theta(frequency)
    #e = 0.0001
    if circ_type =="v0":
        thetas= [CC.gsp_thetas(amps, 3, 9)]
        index=[0,64,128,192,256,320,384]
        thetas.append(CC.indexed_theta(amps,3,9,index))
        index=[0,32,64,128,192,256,320,384]
        for i in range(5):
            thetas.append(CC.indexed_theta(amps,4+i,9,index))
    #e=0.001
    elif circ_type =="v1":
        thetas= [CC.gsp_thetas(amps, 3, 9)]
        index=[0,64,128,256,384]
        for i in range(6):
            thetas.append(CC.indexed_theta(amps,3+i,9,index))
    elif circ_type == "v2":
        thetas= [CC.gsp_thetas(amps, 2, 9)]
        index=[0,128,256]
        for i in range(7):
            thetas.append(CC.indexed_theta(amps,2+i,9,index))
    #e=0.1
    elif circ_type =="v3":
        thetas= [CC.gsp_thetas(amps, 2, 9)]
        index=[0,256]
        for i in range(7):
            thetas.append(CC.indexed_theta(amps,2+i,9,index))
    return thetas


def construct_circ(n,gsp,circ_type,PLOT=True):
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)

    #Reduced frequency range
    frequency = np.linspace(40,170,num=pow(2,n),endpoint=False)
    thetas = gen_angles(frequency,n,circ_type)

    basic = CC.General_State_Prep(thetas[0])
    circ.append(basic,qr[n-gsp:n])

    if circ_type=="v0":
        controls =["000","001","010","011","100","101","11"]
        m_3=CC.ThetaRotation(circ,qr,controls,5,thetas[1],n,True)
        circ.append(m_3,[*qr])
        for i in range((n-4)):
            controls =["0000","0001","001","010","011","100","101","11"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(4+i),thetas[i+2],n,True)
            circ.append(fixed,[*qr])
    elif circ_type=="v1":
        for i in range((n-3)):
            controls =["000","001","01","10","11"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(3+i),thetas[i+1],n,True)
            circ.append(fixed,[*qr])
    elif circ_type=="v2":
        for i in range((n-2)):
            controls =["00","01","1"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(2+i),thetas[i+1],n,True)
            circ.append(fixed,[*qr])
    elif circ_type =="v3":
        for i in range((n-2)):
            controls =["0","1"]
            fixed=CC.ThetaRotation(circ,qr,controls,8-(2+i),thetas[i+1],n,True)
            circ.append(fixed,[*qr])

    circ.save_statevector()
    circ.measure(qr,cla_reg)
    shots = 1000
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    state_vector = result.get_statevector(tqc)
    circ.decompose().draw("mpl",fold=-1)
    plt.show()

    amps = IAFL.amplitudes_for_theta(frequency)
    amps_7_6 = IAFL.amplitudes_7_over_6(frequency)

    #Plots the results in comparision to the actual amplitudes
    if PLOT:
        rc_results = np.sqrt(state_vector.probabilities())
        result_fig = plt.figure()
        result_fig.set_figwidth(8)
        result_fig.set_figheight(6)
        plt.plot(frequency,np.sqrt(state_vector.probabilities()),color='k',label='Statevector',ms = 2)
        plt.plot(frequency,amps_7_6,color='r',label='Amps -7/6',linestyle='dashed')
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        result_fig.savefig("../Images/V0_Inspiral_reduced_frequency.png",dpi=1000,bbox_inches='tight')
        plt.show()

        #Count the number of gates for the circuit
        print(dict(circ.decompose().count_ops()))
        print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))

    print("Fidelity: ",CC.Fidelity(amps_7_6,np.sqrt(state_vector.probabilities())))



if __name__ == "__main__":
    #Correct one
    construct_circ(9,3,"v0")
    #Inspiral_Fixed_Rots(9,3,"v1")
    #Generate thetas
