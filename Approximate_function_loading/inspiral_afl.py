import qiskit as qt
from qiskit.circuit.library.standard_gates import RYGate
import matplotlib.pyplot as plt
import numpy as np
import math
import circ_const as CC

def single_theata(amps,n,m):
    start=0
    increment = int(pow(2,n)/pow(2,m))
    print("m:",m)
    mid =int((0.5)*pow(2,n-m))
    end=int(pow(2,n-m))
    print("mid:",mid)
    print("end:",end)
    upper = sum(amps[start:mid])
    lesser =sum(amps[start:end])
    cos_theta2 = upper/lesser
    cos_theta = np.sqrt(cos_theta2)
    theta = np.arccos(cos_theta)
    return theta*2


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

def calculate_n_k(fmin,fmax,error_rate):
    """
    Args:
    fmin: the minumu frequency of our section
    fmaz: the max frequency of our section
    error_rate: what is the maximum error we are willing to induce

    returns:
    nk: list contaiing eta, and its corresponding k
    """
    #calculate /eta based off Theorem 1 (Marin-Sanchez, et al. 2023)
    lower_deltaf = fmax-fmin
    lower_bel_info = ((0*lower_deltaf)+fmin)**2
    lower_info = np.abs((7/3)*((lower_deltaf**2)/lower_bel_info))
    print("n:",lower_info)

    #Uses eta to calcualte our level (k) where we can switch to fixed bounds
    first_log =np.log(error_rate)
    second_log =math.log((4**-9 - ((96/lower_info**2)*first_log)),2  )
    k = -0.5 * second_log
    print("k:",k)
    nk = [lower_info,k]
    return nk

def fixed_circuit(n,k):
    '''

    :param k:
    :return:
    '''
    k=int(k)
    qr= qt.QuantumRegister(size=n,name='q')
    cla_reg =qt.ClassicalRegister(size=n,name="cla")
    circ = qt.QuantumCircuit(qr,cla_reg)

    frequency = np.linspace(40,170,num=pow(2,n),endpoint=False)
    amps= amplitudes_for_theta(frequency)
    thetas = CC.gsp_thetas(amps,k,n)
    basic=CC.General_State_Prep(thetas)
    circ.append(basic,qr[n-k:n])

    for i in range(n-k):
        theta = single_theata(amps,n,k+i)
        control_y = RYGate(theta)
        circ.append(control_y,qr[(n-(k+i+1)):(n-(k+i))])

    circ.save_statevector()
    #Perform a measurement
    circ.measure(qr,cla_reg)
    shots = 1000
    backend= qt.Aer.get_backend("aer_simulator")
    tqc = qt.transpile(circ,backend)
    job = backend.run(tqc,shots=shots)
    result = job.result()
    state_vector = result.get_statevector(tqc)
    circ.decompose().draw("mpl",fold=-1)
    plt.show()
    print(dict(circ.decompose().count_ops()))
    print(dict(circ.decompose().decompose().decompose().decompose().count_ops()))

    return state_vector

def results(n,state_vector,PLOT=True):
    #Plots the results in comparision to the actual amplitudes
    frequency = np.linspace(40,170,num=pow(2,n),endpoint=False)
    amps= amplitudes_7_over_6(frequency)
    if PLOT:
        rc_results = np.sqrt(state_vector.probabilities())
        result_fig = plt.figure()
        result_fig.set_figwidth(8)
        result_fig.set_figheight(6)
        plt.plot(frequency,np.sqrt(state_vector.probabilities()),color='k',label='Statevector',ms = 2)
        plt.plot(frequency,amps,color='r',label='Goal Inspiral',linestyle='dashed')
        plt.legend()
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
        result_fig.savefig("../Images/Inspiral_AFL_sf_2.png",dpi=1000,bbox_inches='tight')
        plt.show()

        #Count the number of gates for the circuit

        #Plot the box plot to see the errors
        difference = rc_results - amps
        fig = plt.figure()
        fig.set_figwidth(8)
        fig.set_figheight(2)
        #plt.boxplot(abs(difference),vert=False,widths=1)
        plt.scatter(difference,amps)
        #plt.ylim(0,2)
        #plt.yticks([])
        plt.xlabel("Goal")
        plt.ylabel('Real')
        #plt.title("Boxplot showing the errors between the datapoints encoded by the reduced circuit\n and the intended amplitudes.")
        fig.savefig('../Images/Inspiral_AFL_Residuals.png', dpi=1000,bbox_inches='tight')
        plt.show()
    fidelity=CC.Fidelity(amps,np.sqrt(state_vector.probabilities()))
    print("Fidelity: ",fidelity)
    return fidelity


if __name__ == "__main__":
    vals=calculate_n_k(40,170,0.99)

    state_result = fixed_circuit(9,(math.ceil(vals[1])))
    results(9,state_result)