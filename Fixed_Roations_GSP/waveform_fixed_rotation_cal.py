import matplotlib.pyplot as plt
import math
import pickle
import numpy as np



#Provided by Ashwin from LinearSplineApprox.py

def straightline(point1,point2):
    #y=mx+c
    m= (point2[1]-point1[1])/(point2[0]-point1[0])
    c = point1[1]- m*point1[0]
    return [m,c]


def calculate_k(n,error_rate):
    first_log =np.log(error_rate)
    try:
        second_log =np.log2((4**-9 - ((96/n**2)*first_log))  )
    except:
        second_log = 0
    k = -0.5 * second_log
    return k

def rescale_f(current_Freq,fmin,fmax):
    new = (current_Freq - fmin) /(fmax-fmin)
    return new

def numerical_second_devirative(datapoint,fmin,fmax,Plot=True):
    """Get the numerical second log devirative  """
    amps =[]
    x_list=[]
    y_list=[]
    first_devirative=[]
    second_devirative =[]
    count = 1
    for i in datapoint:
        """Rescale the fucntion, then find the log of the function squared"""
        x = rescale_f(i[0],fmin,fmax)
        x_list.append(i[0])
        amps.append(i[1])
        #log of the function squared
        y_list.append(math.log(math.pow(i[1],2)))

    for count,value in enumerate(y_list):
        if count+1!=len(datapoint):
            temp =[x_list[count],value]
            temp_1 = [x_list[count+1],y_list[count+1]]
            stright_line_info =straightline(temp,temp_1)
            #stright_line_info =straightline(i,datapoint[count])
            gradient = stright_line_info[0]
            first_devirative.append(gradient)
        '''else:
            first_devirative.append(gradient)'''

    for count,value in enumerate(first_devirative):
        if count+1 !=len(first_devirative):
            temp =[x_list[count],value]
            temp_1 = [x_list[count+1],first_devirative[count+1]]
            stright_line_info =straightline(temp,temp_1)
            second_devirative.append(stright_line_info[0])

    if Plot:
        fig = plt.figure()
        fig.set_figwidth(8)
        fig.set_figheight(6)
        plt.plot(x_list[:-2], amps[:-2],label="Waveform: f")
        plt.plot(x_list[:-2], y_list[:-2],label="$log(f^2)$",linestyle='dashed')
        plt.plot(x_list[:-2], first_devirative[:-1],label="1st derivative",linestyle='dotted')
        plt.plot(x_list[:-2], second_devirative,label="2nd derivative",linestyle='-.')
        plt.title("Calculating $\eta$ for our numerical distribution")
        plt.xlabel("f (Hz)")
        plt.ylabel("Amp")
        plt.legend()
        fig.savefig('../Images/Waveform_Eta.png', dpi=1000,bbox_inches='tight')
        plt.show()

    eta = max(second_devirative)
    if eta == 0.0:
        eta = second_devirative[0]

    second_devirative.append(stright_line_info[0])
    second_devirative.append(stright_line_info[0])

    return [eta,second_devirative]

'''def numerical_512(n,err,no_datapoint,datafile):
    """
    For 512 datapoints, we calculate eta and corresponding k values for each section
    Note that this code doesn't work for the full waveform as it is not 2nd order differentable.
    Args:
    n: number of qubits/level we want to split up to
    err: the error we are willing to induce
    """
     for 350 datapoints
    with open("datapoints", "rb") as fp:   # Unpickling
        datapoints = pickle.load(fp)
    with open(datafile, 'rb') as pickleFile:
        datapoints = pickle.load(pickleFile)
    fmin =40
    fmax =350
    data_min =0
    no_datapoint = no_datapoint
    print("Testing for err:",err)
    for m in range(n):
        print("level:",m)
        temp_fmin = fmin
        temp_data_min = data_min
        freq_increment = (fmax-fmin)/2**m
        data_increment = no_datapoint/2**m
        for j in range(2**m):
            freq_end = temp_fmin +freq_increment
            data_end = temp_data_min +data_increment
            eta=numerical_second_devirative(datapoints[int(temp_data_min):int(data_end)],temp_fmin,freq_end)[0]
            k=calculate_k(eta,err)
            print("k:",k)
            temp_fmin = freq_end
            temp_data_min = data_end'''


def find_split(split,data):
    """
    Based on the frequency given this return the point in the list which is closest to the split
    :param split: frequency value of the split
    :param data: a list of all the datapoints
    :return: list index
    """
    #print("split",data[0])
    for count,i in enumerate(data):
        if count ==0 and i[0]>split:
            return [None,None]
        elif i[0]>split:
            #print(i[0])
            return [count,i[0]]

    return [None,None]


'''def derivative(data,Plot=True):
    x_list=[]
    y_list=[]
    amps=[]
    for i in data:
       x_list.append(i[0])
       y_list.append(math.log(math.pow(i[1],2)))
       amps.append(i[1])
    y_list =np.array(y_list)
    x=len(y_list)
    div_1=np.gradient(y_list,x)

    div_2 = np.gradient(div_1,x-1)
    if Plot:
        fig = plt.figure()
        fig.set_figwidth(8)
        fig.set_figheight(6)
        plt.plot(x_list[:-2], amps[:-2],label="Waveform: f")
        plt.plot(x_list[:-2], y_list[:-2],label="$log(f^2)$",linestyle='dashed')
        plt.plot(x_list[:-2], div_1[:-2],label="1st derivative",linestyle='dotted')
        #plt.plot(x_list[:-2], div_2,label="2nd derivative",linestyle='-.')
        plt.title("Calculating $\eta$ for our numerical distribution")
        plt.xlabel("f (Hz)")
        plt.ylabel("Amp")
        plt.legend()
        fig.savefig('../Images/Waveform_Eta.png', dpi=1000,bbox_inches='tight')
        plt.show()
    eta = np.abs(div_2)
    return np.max(eta)'''

def divirative_section(current_data,splits,err,fmin,fmax,Plot=True):
    """
    Do the second divitative of 3 different sections, then piece them together
    :param current_data:
    :param splits:
    :return: Divirative of the wholewaveform
    """
    prev_bound =0
    bounds=[]
    results=[]
    x_results=[]
    y_vals=[]
    log_val=[]

    for count,i in enumerate(splits):
        bound = find_split(i,current_data[:])
        if bound[0] == None:
            print("No Bound")
        else:
            bounds.append(bound)
    for count,j in enumerate(current_data[:]):
        x = rescale_f(j[0],fmin,fmax)
        current_data[count][0]=x
        x_results.append(x)
        y_vals.append(j[1])
        log_val.append(np.log(math.pow(j[1],2)))

    for z in bounds:
        temp_data = current_data[prev_bound:z[0]]
        #print("prev:",prev_bound)
        results = results +numerical_second_devirative(temp_data,fmin=fmin,fmax=fmax,Plot=False)[1]
        prev_bound = z[0]
    temp_data= current_data[prev_bound:]
    results = results + numerical_second_devirative(temp_data,fmin=fmin,fmax=fmax,Plot=False)[1]

    if Plot == True:
        plt.plot(x_results, y_vals,label="Waveform")
        plt.plot(x_results, log_val,label="log(f^2(x))",linestyle='-')
        plt.plot(x_results, results,label="2nd derivative",linestyle='-.')
        plt.xlabel("f (Hz)")
        plt.ylabel("Amp")
        plt.legend()
        plt.show()
    abs_result =np.abs(np.array(results))
    eta = max(abs_result)
    print("eta:",eta)
    k=calculate_k(eta,err)
    print("k:",k)
    return results

def waveform_derivative_cal(datafile,n,err):

    fmin =40
    fmax =350
    data_min =0
    no_datapoint = 512
    print("Testing for err:",err)
    splits= [43.728557756818006,243.75275131383347]
    #divirative_section(datapoints,splits,0.99999,40,350,True)

    for m in range(n):
        with open(datafile, 'rb') as pickleFile:
            datapoints = pickle.load(pickleFile)
        print("level:",m)
        temp_fmin = fmin
        temp_data_min = data_min
        freq_increment = (fmax-fmin)/2**m
        data_increment = no_datapoint/2**m
        for j in range(2**m):
            print("j: ",j)
            freq_end = temp_fmin +freq_increment
            data_end = temp_data_min +data_increment
            divirative_section(datapoints[int(temp_data_min):int(data_end)],splits,err=err,fmin= temp_fmin,fmax= freq_end,Plot=True)
            #eta=numerical_second_devirative(datapoints[int(temp_data_min):int(data_end)],temp_fmin,freq_end)[0]

            temp_fmin = freq_end
            temp_data_min = data_end

def cal_eta_k_fixed_bound(datafile,n,error_rate):
    fmin =40
    fmax =350
    print("Testing for err:",error_rate)
    splits= [43.728557756818006,243.75275131383347]
    with open(datafile, 'rb') as pickleFile:
        datapoints = pickle.load(pickleFile)
    full_eta_result = np.array( divirative_section(datapoints,splits,err=error_rate,fmin= fmin,fmax= fmax,Plot=True))
    print("\n----------\n")
    for m in range(n):
        print("level:",m)
        temp_fmin = int(0)
        increment = int((512-0)/2**m)
        for j in range(2**m):
            end = temp_fmin +increment
            #Here add the factor
            eta = math.pow(2,-(m)) *full_eta_result[temp_fmin:end]
            max_eta = max(np.abs(eta))
            print("\nn:",max_eta)
            k=calculate_k(max_eta,error_rate)
            print("k:",k)
            temp_fmin = end
        print("\n----------\n")



if __name__ == "__main__":
    #insprial_data =inspiral_point()
    datafile="waveform_data.pkl"
    #divirative_section(datafile,[43.728557756818006,243.75275131383347],0.99999)
    cal_eta_k_fixed_bound(datafile,5,0.9999)
    #waveform_derivative_cal(datafile,5,0.9999)
    #numerical_512(3,0.9999,512,datafile)