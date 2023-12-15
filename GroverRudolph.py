'''This program is disigned to run a generic version of the Grover-Rudolph algorithm
It should take a probabiltiy distribution and encoded it into a quantum state of m qubits'''
import qiskit as qt
import matplotlib.pyplot as plt
import numpy as np
import scipy

def prob(qubit,fmin,fdelta,distrbution_type):
    '''
    Args:
    qubit: takes in the Qubit (we will use this in a decimal form)
    fmin: minimum frequency we are using
    fdelta: the change in frequence per each step
    distribution_type: the number of the power to be applied to this
    Note in Hayes paper there distribution_type is 7/3
    Return:
    probabilty from that point in the qubit :p(i) in the 
    '''
    '''f=fmin + i delta f
    so p(i) = (fmin+i delta f)^(=7/3) / Nc
    Where Nc is the normalisation constant
    delta f = (fmin-fmax)/2^n '''
    #print(qubit)
    probabilityI = (fmin + (qubit*fdelta))**distrbution_type
    #This is where we need the probability distribution
    print("prob",probabilityI)
    return probabilityI

def cos2_theta(m,j,n,distribution):
    '''Args:
    m: the level we are at
    j: which bin we are in
    n:
    distribution:
    Return:
    The angle theta 
    '''
    #two sums which are then divided over each other
    i = j*2**(n-m-1)
    #I think this is where i need to minus 1 from the equation
    #Think about it for the sum though?
    half = (j+0.5)*2**(n-m-1)
    full =  (j+1)*2**(n-m-1)
    #print(i)
    #Calculate the constants from the classical distribiution for prob
    fmin = distribution[0]
    size = len(distribution)
    fmax = distribution[(size-1)]
    #normalise our distribution
    normal = np.sqrt(np.sum(np.abs(distribution)**2))
    print("normal:",normal)
    fdelta = (fmax-fmin)/(2**n)
    print("fdelta",fdelta)
    top =0
    bottom =0
    #CHANGE FOR THE DIFFERENT KINDS OF PROBABILITY
    #THIS CURRENTLY DOES NOTHING
    distrbution_type =1
    for x in range(i,int(half)):
        #print(x)
        top = top + prob(x,fmin,fdelta,distrbution_type)
    
    for y in range(i,full):
        bottom = bottom + prob(y,fmin,fdelta,distrbution_type)
    
    #print(top)
    #print(bottom)
    
    result = top/bottom

    #Now get the angle from the cos^2(theta)

    return result

def Grover_Rudolph_func(n,distribution):
    '''Args: 
    n: the number of qubits for our circuit
    distribution: is the probability distrbution we are trying to encode
    Returns: an encoded quantum circuit '''

    #Define how many levels we want to have
    m = list(range(0, n-1))
    #m=[0,1]
    #Initalise the quantum circuit
    qc = qt.QuantumCircuit(n)
    angles={}
    #Loop through each level of m
    for i in m:
        #split up the probability distribution into two parts
        #print("m",i)
        #This is our j in the maths
        #Define the number of bins j (how many bins we split the probability distribution into)
        current_bins = list(range(0,2**i))
        #print(current_bins)
        for j in current_bins:
            #print("j",j)
            theta=cos2_theta(i,j,n,distribution)
            place=str(i)+str(j)
            angles[place]=theta
    print(angles)
    #Draws the quantum circuit
    #qc.draw("mpl")
    #plt.show()

if __name__ == "__main__":
    Grover_Rudolph_func(3,[0,1,2,3,4,5,6,7])
