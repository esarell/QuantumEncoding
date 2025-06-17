'''Just putting all my code for General State Prepration in One place.
This is the Grover-Rudolph Method for amplitude encoding'''
import pickle
import numpy as np
from numpy.ma.core import absolute
from waveform_fixed_rotation_cal import calculate_k

def waveform_amps():
    """

    :return:
    """
    with open("waveform_amps.pkl", 'rb') as pickleFile:
        test = pickle.load(pickleFile)
    return test

def prob_normalised(data):
    """

    :param data:
    :return:
    """
    result=[]
    for x in data:
        square = x*x
        norm = np.sqrt(square/sum(data*data))
        result.append(norm)
    return result

def waveform_theta(frequency,m):
    """
    :param frequency:
    :param m:
    :return:
    """
    j=0
    #This is specifically for f^-7/6 (so p(i)is f^-7/3)
    #Will need changing for a more complicatated distribution
    #Gets a list of normalised amplitudes for all the frequencies
    temp_amps= waveform_amps()
    amps=prob_normalised(temp_amps)
    thetas=[]
    for i in range(m):
        j=pow(2,i)
        gsp_theta=[]
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
        thetas.append(gsp_theta)
    return thetas


def sig_square(k):
    return (1/pow(2,(k-1)))

def Compare_Theata(n,Equation=17,RC=False):
    '''
    Compares the difference between the two theta values to see the eta value associated with it.
    Based on equation B17 by Marin-Sanchez et al.
    :param n:
    :return:
    '''
    n=9
    frequency = np.linspace(40,350,num=pow(2,n),endpoint=False)
    angles = waveform_theta(frequency,n)
    if Equation == 21:
        print("Equation B21")
        for count,i in enumerate(angles):
            if count ==0:
                print("no comparison for level m=0")
                '''sig = 1/pow(2,1)
                first = (abs([2]-i[1])*4)/pow(sig,2)
                print("m 1:",first)'''
            else:
                #Level 2
                sig = sig_square(count)
                absolute_diff = (abs(i[0]-i[-1]))
                eta_approx = (absolute_diff/sig)*4
                print("m ",count)
                print("diff",absolute_diff)
                print("eta",eta_approx)
                print ("k",calculate_k(eta_approx,0.9999))
                print("\n-----------\n")
    elif Equation == 22:
        print("Equation B22")
        for count,i in enumerate(angles):
            #This is the angle differences for my circuit
            if RC==True:
                if count ==0:
                    print("no comparison for level m=0")
                elif count == 1 or count ==2:
                    sig = sig_square(count)
                    total =0
                    for j,x in enumerate(i):
                        try:
                            absolute_diff = (abs(x-i[j+1]))
                            total = total + absolute_diff
                        except:
                            total = total
                    eta_approx = (total/sig)*4
                    print("m ",count)
                    print("total_diff",absolute_diff)
                    print("eta",eta_approx)
                    print ("k",calculate_k(eta_approx,0.9999))
                    print("\n-----------\n")
                elif count == 3:
                    total=0
                    sig = sig_square(count)
                    diff_results=[]
                    eta_result=[]
                    k_result=[]
                    print("m ",count)
                    splits = ["000","001","010","011","10","11"]
                    print(len(i[:-4]))
                    for j,x in enumerate(i[:-1]):
                        absolute_diff = (abs(x-i[j+1]))
                        total = total + absolute_diff
                        if j >= 4:
                            diff_results.append(absolute_diff)
                            eta_approx = (total/sig)*4
                            eta_result.append(eta_approx)
                            k_result.append(calculate_k(eta_approx,0.9999))
                            total=0
                    print("m ",count)
                    print("total_diff",diff_results)
                    print("eta",eta_result)
                    print ("k",k_result)
                    print("\n-----------\n")
                    absolute_diff=0
                elif count >= 4:
                    print("m ",count)
                    splits = [0,0.0625,0.125,0.25,0.375,0.5,0.75,1]
                    for s_count,s in enumerate(splits[:-1]):
                        total=0
                        start = int(s*len(i))
                        end = int(splits[s_count+1]*len(i))
                        for j,x in enumerate(i[start:end]):
                            absolute_diff = (abs(x-i[j+1]))
                            total = total + absolute_diff
                        eta_approx = (total/sig)*4
                        print("total_diff",absolute_diff)
                        print("eta",eta_approx)
                        print ("k",calculate_k(eta_approx,0.9999))
                    print("\n-----------\n")


            else:
                if count ==0:
                    print("no comparison for level m=0")
                    '''sig = 1/pow(2,1)
                    first = (abs([2]-i[1])*4)/pow(sig,2)
                    print("m 1:",first)'''
                else:
                    #Level 2
                    sig = sig_square(count)
                    total =0
                    for j,x in enumerate(i):
                        try:
                            absolute_diff = (abs(x-i[j+1]))
                            total = total + absolute_diff
                        except:
                            total = total
                    eta_approx = (total/sig)*4
                    print("m ",count)
                    print("total_diff",absolute_diff)
                    print("eta",eta_approx)
                    print ("k",calculate_k(eta_approx,0.9999))
                    print("\n-----------\n")



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

if __name__ == "__main__":
    Compare_Theata(9,22)
    Compare_Theata(9,22,RC=True)