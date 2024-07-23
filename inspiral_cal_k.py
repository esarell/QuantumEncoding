import matplotlib.pyplot as plt
import math

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
    error_rate =0.999
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
    print("error_rate:",error_rate)
    lower_deltaf = fsplit-fmin
    lower_bel_info = ((0*lower_deltaf)+fmin)**2
    lower_info = (7/3)*((lower_deltaf**2)/lower_bel_info)
    print("n:",lower_info)
    first_log =math.log(error_rate,10)
    second_log =math.log((4**-9 - ((96/lower_info**2)*first_log)),2  )
    k = -0.5 * second_log
    print("k:",k)





#result=optimal_fsplit()
#cal_k(result[0])
#cal_k(result[1])
print("test 1:no controls")
calculate_n_k(40,350,0.9999)
print("\ntest 2: 1 controls")
calculate_n_k(40,195,0.9999)
calculate_n_k(195,350,0.9999)
#calculate_n_k(155,253,0.9997)
#calculate_n_k(98,155,0.9997)
#calculate_n_k(253,350,0.9997)
print("\ntest 3: 2 controls")
calculate_n_k(40,117.5,0.9999)
calculate_n_k(117.5,195,0.9999)
calculate_n_k(195,272.5,0.9999)
calculate_n_k(272.5,350,0.9999)

print("\ntest 4: 3 controls")
calculate_n_k(40,78.75,0.9999)
calculate_n_k(78.75,117.5,0.9999)
calculate_n_k(117.5,156.25,0.9999)
calculate_n_k(156.25,195,0.9999)
calculate_n_k(195,233.75,0.9999)
calculate_n_k(233.75,272.5,0.9999)
calculate_n_k(272.5,311.25,0.9999)
calculate_n_k(311.25,350,0.9999)

print("\ntest 5: 4 controls\njust for first bits")
calculate_n_k(40,59.375,0.9999)
calculate_n_k(59.375,78.75,0.9999)
calculate_n_k(78.75,98.125,0.9999)
calculate_n_k(98.125,117.5,0.9999)

print("Conclusion: dont need more than 4 controls")
#calculate_n_k(40,69,0.99999)
#calculate_n_k(69,98,0.99999)
#calculate_n_k(98,127,0.99999)
#calculate_n_k(127,155,0.99999)
#calculate_n_k(155,204,0.99999)
#calculate_n_k(204,253,0.99999)
#calculate_n_k(253,302,0.99999)
#calculate_n_k(302,350,0.99999)

#calculate_n_k(40,55,0.9999999)
#calculate_n_k(55,69,0.9999999)
#calculate_n_k(69,84,0.9999999)
#calculate_n_k(84,98,0.9999999)
#calculate_n_k(98,113,0.9999999)
#calculate_n_k(113,127,0.9999999)
#calculate_n_k(127,141,0.9999999)
#calculate_n_k(141,155,0.9999999)




#result_n=plot_n()
#cal_k(result_n)
#result_fmin=plot_n_movingfmin()
#cal_k(result_fmin)