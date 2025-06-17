import numpy as np
import matplotlib.pyplot as plt
#Code for the Black-Scholes distribution based on Eq. 35 from Marin-Sachez et al. 2023

def calculate_k(n,eta,error_rate):
    first_log =np.log10(error_rate)
    try:
        second_log =np.log2((4**-n - ((96/eta**2)*first_log))  )
    except:
        second_log = 0
    k = -0.5 * second_log
    return k

def rescale(current_x,x_min,x_max):
    new = (current_x - x_min) /(x_max-x_min)
    return new

def straightline(point1,point2):
    #y=mx+c
    m= (point2[1]-point1[1])/(point2[0]-point1[0])
    c = point1[1]- m*point1[0]
    return [m,c]

def factor_eta_cal():
    #We change the bounds and then apply a factor of 2^-m
    return 0

def eta_cal(x,y,x_min,x_max,Plot=True):
    """Get the numerical second log devirative  """
    amps =[]
    x_list=[]
    y_list=[]
    first_devirative=[]
    second_devirative =[]
    count = 1
    for count,i in enumerate(x):
        """Rescale the fucntion, then find the log of the function squared"""
        x_list.append( rescale(i,x_min,x_max))
        #log of the function squared
        y_list.append(np.log(np.square(y[count])))

    for count,value in enumerate(y_list):
        if count+1!=len(x):
            temp =[x_list[count],value]
            temp_1 = [x_list[count+1],y_list[count+1]]
            stright_line_info =straightline(temp,temp_1)
            gradient = stright_line_info[0]
            first_devirative.append(gradient)

    for count,value in enumerate(first_devirative):
        if count+1 !=len(first_devirative):
            temp =[x_list[count],value]
            temp_1 = [x_list[count+1],first_devirative[count+1]]
            stright_line_info =straightline(temp,temp_1)
            second_devirative.append(np.abs(stright_line_info[0]))

    if Plot:
        print(len(x_list))
        print("amps",len(y_list))
        fig = plt.figure()
        fig.set_figwidth(8)
        fig.set_figheight(6)
        plt.plot(x_list[:], y[:],label="BS")
        plt.plot(x_list[:], y_list[:],label="$log(f(x)^2)$",linestyle='dashed')
        plt.plot(x_list[:-1], first_devirative[:],label="1st derivative",linestyle='dotted')
        plt.plot(x_list[:-2], second_devirative,label="2nd derivative",linestyle='-.')
        plt.title("Calculating $\eta$ for our numerical distribution")
        plt.xlabel("x")
        plt.ylabel("Amplitudes")
        plt.legend()
        fig.savefig('../Images/BS_Eta.png', dpi=1000,bbox_inches='tight')
        plt.show()

    eta = max(second_devirative)
    if eta == 0.0:
        eta = second_devirative[0]
    second_devirative.append(stright_line_info[0])
    second_devirative.append(stright_line_info[0])

    return [eta,second_devirative]

def rescale_eta(n,x,y,nCal=6):
    x_min = -10
    x_max = 10
    for m in range(nCal):
        print("level:",m)
        current_x = x_min
        start_place = 0
        x_increment = (x_max-x_min)/2**m
        y_increment = len(y)/2**m
        for j in range(2**m):
            print("j: ",j)
            x_end = current_x +x_increment
            end_place = start_place +y_increment
            eta =eta_cal(x[int(start_place):int(end_place)],y[int(start_place):int(end_place)],current_x,x_end,Plot=False)[0]
            #eta=numerical_second_devirative(datapoints[int(temp_data_min):int(data_end)],temp_fmin,freq_end)[0]
            print("eta:",eta)
            print("k:",calculate_k(n,np.abs(eta),0.9999))
            current_x = x_end
            start_place = end_place


def Black_Scholes(x,K,s):
    #Genrates the distribution
    y = []
    test = np.log(K*s)
    print(-test)
    for i in x:
        if -test <= i < 0:
            y.append( K-(np.exp(-i)/s))
        elif 0 < i <= test:
            y.append( K-(np.exp(i)/s))
        else:
            y.append(0.0000001)
    return y

def Gen_BS(n,Plot=True):
    K = 45
    c = 3
    s = K*c
    x_values = np.linspace(-8,8,num=pow(2,n),endpoint=True)
    y_values = Black_Scholes(x_values,K,s)
    if Plot:
        plt.plot(x_values,y_values)
        plt.xlabel("x")
        plt.ylabel("Amplitudes")
        plt.show()
    rescale_eta(n,x_values,y_values)



if __name__ == "__main__":
    Gen_BS(8,Plot=False)

