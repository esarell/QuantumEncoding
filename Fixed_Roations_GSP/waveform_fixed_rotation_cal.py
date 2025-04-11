import matplotlib.pyplot as plt
import math
import pickle
import csv
import numpy as np
from mpmath.libmp.libelefun import log_taylor_cache


#Provided by ashwin from LinearSplineApprox.py 33 points
#datapoints = [[40.0, 45.53588379231391], [41.73988252356714, 42.99468907854167], [43.60341387652148, 40.534160560597044], [45.60224410550622, 38.15604254353809], [47.750013001775706, 35.86154967907757], [50.06179757799808, 33.65242668859051], [52.55504724854148, 31.530155680572275], [55.2530432923146, 29.49343653352211], [58.18193853600308, 27.541866623249554], [61.36855354444498, 25.677762482084642], [64.8414792064986, 23.90506569350602], [68.63900132055133, 22.225094786058857], [72.80994794058644, 20.637790275494588], [77.4080959125527, 19.145254748772437], [82.50040406939527, 17.74900736770979], [88.17265962473014, 16.449716530123492], [94.52432260915532, 15.25000778588201], [101.68458619910827, 14.15186546357962], [109.82595132795439, 13.15617543160781], [119.21421437293974, 12.259090629726577], [130.20746393759126, 11.457567247773428], [143.36827004135537, 10.743344958945583], [159.8044251378891, 10.090026336348101], [182.26246814673658, 9.419528957212911], [225.86535582759427, 8.224181668209477], [259.0729889700747, 6.826576178869628], [269.37356128061754, 5.967553107437675], [280.82588990553296, 4.772774132870048], [309.22430385343546, 2.0406806963843023], [320.6491715928424, 1.357629093553048], [333.2701260238657, 0.8568326977345254], [348.555883309873, 0.4927565704511884], [350.0, 0.46802388200131884]]

def straightline(point1,point2):
    #y=mx+c
    m= (point2[1]-point1[1])/(point2[0]-point1[0])
    c = point1[1]- m*point1[0]
    return [m,c]

def gen_lin_vals(linearParams, fmin,fmax,no_point):
    increment = (fmax-fmin)/no_point
    data=[]
    for i in no_point:
        y=linearParams[0]*fmin+linearParams[1]
        data.append([fmin,y])
        fmin= fmin+increment
    return data

def inspiral_sec_log_der(a,x,b,fmin,fmax):
    """Args:
    a:the gradient of the straight line
    b:the y intercept from the staright line
    x: value we are evaluating
    fmin:the starting frequency of the line
    fmax:the ending freqency of the line  """
    denominator = math.pow(((a*((x*(fmax-fmin))+fmin))+b),2)
    n = (-2*math.pow(a,2))/denominator
    return n


def inspiral_fun(f):
    amp = math.pow(f,(-7/6))
    return amp

def inspiral_point():
    no =33
    fmin =40
    fmax =350
    delta = (fmax-fmin)/no
    data=[]
    for i in range(34):
        amplitude =inspiral_fun(fmin)
        result=[fmin,amplitude]
        data.append(result)
        fmin=fmin+delta
    return data

def cal_n_1_over_log(x,fmin,fmax):
    delta = fmax-fmin
    top_frac = 2*math.pow(delta,2)*(math.log((delta*x)+fmin)+1)
    bottom_frac = math.pow(delta*x+fmin,2)*math.pow(math.log(delta*x+fmin),2)
    n = top_frac/bottom_frac
    return n 

def cal_anltical_inspiral(x,fmax,fmin):
    deltaf = fmax-fmin
    bel_info = ((x*deltaf)+fmin)**2
    info = (7/3)*((deltaf**2)/bel_info)

def calculate_k(n,error_rate):
    first_log =math.log(error_rate,10)
    try:
        second_log =math.log((4**-9 - ((96/n**2)*first_log)),2  )
    except:
        second_log = 0
    k = -0.5 * second_log
    return k

def rescale_f(current_Freq,fmin,fmax):
    new = (current_Freq - fmin) /(fmax-fmin)
    return new



def numerical_second_devirative(datapoint,Plot=False):
    """Get the numerical second log devirative  """
    amps =[]
    x_list=[]
    y_list=[]
    first_devirative=[]
    second_devirative =[]
    count = 1
    for i in datapoint:
        #x = rescale_f(i[0],fmin,fmax)
        x_list.append(i[0])
        amps.append(i[1])
        #log of the function squared
        y_list.append(math.log(math.pow(i[1],2)))

    print(len(y_list))
    for count,value in enumerate(y_list):
        if count+1!=len(datapoint):
            temp =[x_list[count],value]
            temp_1 = [x_list[count+1],y_list[count+1]]
            stright_line_info =straightline(temp,temp_1)
            '''stright_line_info =straightline(i,datapoint[count])
            gradient = stright_line_info[0]
            first_devirative.append(gradient)'''
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
    print("eta:",eta)
    print("len:",len(second_devirative))
    second_devirative.append(stright_line_info[0])
    second_devirative.append(stright_line_info[0])
    return second_devirative

def numerical_512(n,err,no_datapoint,datafile):
    """
    For 512 datapoints, we calculate eta and corresponding k values for each section 
    Args:
    n: number of qubits/level we want to split up to
    err: the error we are willing to induce
    """
    ''' for 350 datapoints
    with open("datapoints", "rb") as fp:   # Unpickling
        datapoints = pickle.load(fp)'''
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
            eta=numerical_second_devirative(datapoints[int(temp_data_min):int(data_end)],temp_fmin,freq_end,Plot=True)
            k=calculate_k(eta,err)
            print("k:",k)
            temp_fmin = freq_end
            temp_data_min = data_end

def find_split(split,data):
    """

    :param split:
    :param datapoints:
    :return: list index
    """
    print("split",data[0])
    for count,i in enumerate(data):
        if i[0]>split:
            print(i[0])
            return [count,i[0]]


def derivative(data,Plot=True):
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
    return np.max(eta)

def divirative_section(datafile,splits,err):
    """
    Do the second divitative of 3 different sections, then piece them together
    :param datafile:
    :param splits:
    :return: Divirative of the wholewaveform
    """
    prev_bound =0
    fmin =40
    with open(datafile, 'rb') as pickleFile:
        datapoints = pickle.load(pickleFile)
    print(datapoints[0])
    bounds=[]
    results=[]
    x_results=[]
    y_vals=[]
    log_val=[]
    '''eta = derivative(datapoints,Plot=True)
    print("eta:",eta)
    k=calculate_k(eta,err)
    print("k:",k)'''
    for count,i in enumerate(splits):
        bound = find_split(i,datapoints)
        bounds.append(bound)

    for count,j in enumerate(datapoints):
        x = rescale_f(j[0],40,350)
        datapoints[count][0]=x
        x_results.append(x)
        y_vals.append(j[1])
        log_val.append(math.log(math.pow(j[1],2)))

    for z in bounds:
        temp_data = datapoints[prev_bound:z[0]]
        print("prev:",prev_bound)
        results = results +numerical_second_devirative(temp_data,Plot=True)
        prev_bound = z[0]
    temp_data= datapoints[prev_bound:]
    results = results + numerical_second_devirative(temp_data,Plot=True)

    print("here")
    plt.plot(x_results, y_vals,label="Waveform")
    plt.plot(x_results, log_val,label="log(f^2(x))",linestyle='-')
    plt.plot(x_results, results,label="2nd derivative",linestyle='-.')
    plt.xlabel("f (Hz)")
    plt.ylabel("Amp")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    #insprial_data =inspiral_point()
    datafile="waveform_data.pkl"
    divirative_section(datafile,[43.728557756818006,243.75275131383347],0.99999)
    #numerical_512(3,0.9999,512,datafile)