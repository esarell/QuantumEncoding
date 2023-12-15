'''
This file is just going to be used as a very simple was of generating different
probability distributions to test the GroverRudolph algorithm
There will be a straight line,
exponentail,
gaussian.
Hayes tested on f^(-7/3) with fmin=40Hz, fmax=168Hz, fdelta = 2Hz, requiring 6 qubits Tdelta =0.02s
'''
def straightLine(n):
    '''Args:
        n is the number of quibits the system is using
        Returns:
        A list of N=2^n amplitudes, creating a probability distribution of (x=y)
    '''
    #This is not currently normalised
    amplitudes=[]
    for i in range((2**n)):
        amplitudes.append(i)
    print(amplitudes)
    return amplitudes

if __name__ == "__main__":
    print("Please use this file as a libary for the GroverRudoph.py file")
