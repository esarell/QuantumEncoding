'''This program is disigned to run a generic version of the Grover-Rudolph algorithm
It should take a probabiltiy distribution and encoded it into a quantum state of m qubits'''
import qiskit as qt
import matplotlib.pyplot as plt
import numpy as np

def cos2_theta (m,j):
    '''Args:
    m: the level we are at
    j: which bin we are in
    Return:
    The angel theta 
    '''


def Grover_Rudolph_func(n,prob):
    '''Args: 
    n: the number of qubits for our circuit
    prob: is the probability distrbution we are trying to encode
    Returns: an encoded quantum circuit '''

    #Define how many levels we want to have
    m = list(range(0, n-1))
    #Initalise the quantum circuit
    qc = qt.QuantumCircuit(n)

    #Loop through each level of m
    for i in m:
        #split up the probability distribution into two parts
        print(i)
        #This is our j in the maths
        #Define the number of bins j (how many bins we split the probability distribution into)
        current_bins = list(range(0,2**i))
        print(current_bins)

    #Draws the quantum circuit
    qc.draw("mpl")
    plt.show()

if __name__ == "__main__":
    Grover_Rudolph_func(8,3)
