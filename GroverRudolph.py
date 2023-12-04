'''This program is disigned to run a generic version of the Grover-Rudolph algorithm
It should take a probabiltiy distribution and encoded it into a quantum state of m qubits'''
import qiskit as qt
import matplotlib.pyplot as plt
import numpy as np
#test
def prob(qubit):
    '''
    Args:
    takes in the place? Qubit?
    Return:
    The probabibility of the qubit
    '''
    #prob qubtit(0,0)
    #|\bra{state}\ket{a}|^2
    #Do we have a valule a or is just probabilty of p000

    return qubit


def cos2_theta(m,j,n):
    '''Args:
    m: the level we are at
    j: which bin we are in
    Return:
    The angel theta 
    '''
    #two sums which are then divided over each other
    i = j*2**(n-m-1)
    #I think this is where i need to minus 1 from the equation
    #Think about it for the sum though?
    half = (j+0.5)*2**(n-m-1)
    full =  (j+1)*2**(n-m-1)
    #print(i)
    top =0
    bottom =0
    for x in range(i,int(half)):
        #print(x)
        top = top + prob(x)
    
    for y in range(i,full):
        bottom = bottom + prob(y)
    
    #print(top)
    #print(bottom)
    
    result = top/bottom

    #Now get the angle from the cos^2(theta)

    return result

def Grover_Rudolph_func(n,prob):
    '''Args: 
    n: the number of qubits for our circuit
    prob: is the probability distrbution we are trying to encode
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
            theta=cos2_theta(i,j,n)
            place=str(i)+str(j)
            angles[place]=theta
    print(angles)
    #Draws the quantum circuit
    #qc.draw("mpl")
    #plt.show()

if __name__ == "__main__":
    Grover_Rudolph_func(3,3)
