import matplotlib.pyplot as plt
import math
import numpy as np
'''This file calculates the relevant \eta values and corresponding k values
for the inspiral of our gravitational wave approx f^{7/6} '''

def plot_n():
    value = 0
    fmin =40
    #fmax =350
    xpoints =[]
    ypoints =[]
    for fmax in range (41,350):
        deltaf = fmax-fmin
        bel_info = ((0*deltaf)+fmin)**2
        info = (7/3)*((deltaf**2)/bel_info)
        print(fmax,info)
        delta = 1
        value = value + delta
        xpoints.append(fmax)
        ypoints.append(info)
    plt.axhline(0,color='red')
    plt.axhline(8*math.pi,color='blue')
    plt.plot(xpoints, ypoints)
    plt.xlabel("fmax")
    plt.ylabel("n")
    plt.show()
    data = [xpoints,ypoints]
    return data    

def plot_n_movingfmin():
    value = 0
    fmin =40
    fmax =350
    xpoints =[]
    ypoints =[]
    for fmin in range (40,349):
        deltaf = fmax-fmin
        bel_info = ((0*deltaf)+fmin)**2
        info = (7/3)*((deltaf**2)/bel_info)
        print(fmin,info)
        delta = 1
        value = value + delta
        xpoints.append(fmin)
        ypoints.append(info)
    plt.axhline(20,color='red')
    plt.axhline(0,color='red')
    plt.plot(xpoints, ypoints)
    plt.xlabel("fmin")
    plt.ylabel("n")
    plt.show()
    data = [xpoints,ypoints]
    return data

def cal_k(data):
    xpoints =data[0]
    ypoints =data[1]
    error_rate =0.9999
    k_values =[]
    first_log =math.log(error_rate,10)
    for n in ypoints:
        second_log =math.log((4**-9 - ((96/n**2)*first_log)),2  )
        k = -0.5 * second_log
        print(n,k)
        k_values.append(k)
    plt.axvline(9,color='red')
    plt.axvline(2,color='red')
    plt.axhline(8*math.pi,color='blue')
    plt.plot(k_values,ypoints)
    plt.xlabel("k")
    plt.ylabel("n")
    plt.show()   

def optimal_fsplit():
    """
    Returns: a 2 dimensional array containg two sets of coordinates
    one for the lower point one for our upper points
    """
    lower_bound = 91
    upper_bound = 171
    fmin =40
    fmax = 350
    lower_xpoints =[]
    lower_ypoints =[]

    upper_xpoints =[]
    upper_ypoints =[]
    for i in range(91,171):
        lower_deltaf = i-fmin
        #why is x =0 double check 
        #I think this was giving the maximum
        lower_bel_info = ((0*lower_deltaf)+fmin)**2
        lower_info = (7/3)*((lower_deltaf**2)/lower_bel_info)
        print(fmin,lower_info)
        lower_xpoints.append(i)
        lower_ypoints.append(lower_info) 
        upper_deltaf = fmax-i
        upper_bel_info = ((0*upper_deltaf)+i)**2
        upper_info = (7/3)*((upper_deltaf**2)/upper_bel_info)
        upper_xpoints.append(i)
        upper_ypoints.append(upper_info) 

    plt.plot(lower_xpoints, lower_ypoints)
    plt.plot(upper_xpoints, upper_ypoints)
    plt.axhline(0,color='red')
    plt.axhline(8*math.pi,color='blue')
    plt.legend(['fmin =40','fmax = 350'])
    plt.xlabel("fsplit")
    plt.ylabel("n")
    plt.show()
    data = [[lower_xpoints,lower_ypoints],[upper_xpoints,upper_ypoints]]
    return data

def label_gate_concept():
    split_value = 118
    fmin =[40]
    fmax =[350]
    half = (fmax-fmin)/2
    fmin.append(half)
    fmax.append(half)

    quater = (fmax[0]-fmin[1])/2
    quater2 = (fmax[1] - fmin[0])/2
    

def calculate_n_k(fmin,fsplit,error_rate):
    """
    Args:
    fmin: the minumu frequency of our section
    fsplit: the max frequency of our section
    error_rate: what is the maximum error we are willing to induce

    returns:
    nk: list contaiing eta, and its corresponding k
    """
    #calculate /eta based off Theorem 1 (Marin-Sanchez, et al. 2023)
    lower_deltaf = fsplit-fmin
    lower_bel_info = ((0*lower_deltaf)+fmin)**2
    lower_info = (7/3)*((lower_deltaf**2)/lower_bel_info)
    print("n:",lower_info)

    #Uses eta to calcualte our level (k) where we can switch to fixed bounds
    first_log =math.log(error_rate,10)
    second_log =math.log((4**-9 - ((96/lower_info**2)*first_log)),2  )
    k = -0.5 * second_log
    print("k:",k)
    nk = [lower_info,k]
    return nk

def cal_eta_k_fixed_bound(n,error_rate):
    '''
    Instead of rescaling we are trying to add a factor see what the outcome is
    :param n:
    :return:
    '''
    values = np.linspace(0,1,num=pow(2,9),endpoint=True)
    print(len(values))
    lower_deltaf = 350-40
    lower_bel_info = ((values*lower_deltaf)+40)**2

    lower_info = (7/3)*((lower_deltaf**2)/lower_bel_info)
    for m in range(n):
        print("level:",m)
        temp_fmin = int(0)
        increment = int((512-0)/2**m)
        print("inc:",increment)
        print(len(lower_info))
        for j in range(2**m):
            end = temp_fmin +increment
            #Here add the factor
            eta = math.pow(2,-(m)) *lower_info[temp_fmin]
            print("\nn:",eta)
            first_log =math.log(error_rate,10)
            second_log =math.log((4**-9 - ((96/eta**2)*first_log)),2  )
            k = -0.5 * second_log
            print("k:",k)
            temp_fmin = end
        print("\n------\n")




    #Uses eta to calcualte our level (k) where we can switch to fixed bounds





if __name__ == "__main__":
    #Splits up my distribution for n levels 
    #Calculates a eta and corresponding k value
    err = 0.9999
    '''n = 3
    fmin =40
    fmax = 350
    print("Testing for err:",err)
    for m in range(n):
        print("level:",m)
        temp_fmin = fmin
        increment = (fmax-fmin)/2**m
        for j in range(2**m):
            end = temp_fmin +increment
            calculate_n_k(temp_fmin,end,err)
            temp_fmin = end'''
    cal_eta_k_fixed_bound(6,err)