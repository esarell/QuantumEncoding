import numpy as np
import matplotlib.pyplot as plt
import scipy

def g(x1): #This is the function
    return np.arccos(x1)

def g2(x): #This is the second derivative of g(x)
    d2 = -x/(1-(x**2))**(3/2)
    return d2

def gen_bounds(n): 
    '''Generates the bounds depended on how many qubits n we have
    args: n: number of qubits
    returns: '''
    k = 4.5e-3  # n = 4
    #k = 1.7e-3  # n = 5
    #k = 0.48e-3 # n = 6 
    if n ==4:
       k = 4.5e-3
    elif n ==5:
        k = 1.7e-3  
    #points are the x axis
    #g(points) is y axis
    points = np.array([])
    xi = 0.01
    xn = xi
    points = np.append(points,xi)
    while xn < 0.99:
        delx = (k/np.abs(g2(xn)))**(1/2)
        xn = xn + delx
        points = np.append(points,xn)
    
    data=[0,1]
    data[0]=points
    data[1]=g(points)    
    print(len(points))
    print(data)

    return data
