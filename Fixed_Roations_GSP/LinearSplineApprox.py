#!/usr/bin/env python
# coding: utf-8

# In[1]:

import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, argparse
import scipy


# In[2]:


matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
width=0.75
color='black'
fontsize=28
ticksize=22
figsize=(10,10)


# # Fergus's work

# In[3]:


phase = True
fmin=40.
fmax=168.
m1=35.
m2=30.
Tfrac=100.
beta=0.
sig=0.
nx = 6
nintx = int(np.ceil(np.log2((fmax-fmin))))

xmax = np.power(2,nintx) - np.power(2,nintx-nx)
xmin = 0.

def m_geo(m):
    return (4.926e-6)*m

df = (fmax-fmin)/(2**nx)
T = 1./df

####### Physical system parameters ###################
m1 = m_geo(m1)
m2 = m_geo(m2)
tc = T + (T/Tfrac)
DT = tc%T
Mt = m1 + m2
nu = (m1*m2)/Mt
eta = nu/Mt
Mc = Mt*eta**(3./5)

def x_trans(x):
    x = x/xmax
    x = x*(fmax-fmin-df)
    x = x + fmin
    return x

xs = np.linspace(xmin,xmax,2**(nx))
xsT = x_trans(xs)

def f_x(x, eps1=1, eps2=1, factor=1):
    x = x_trans(x)
    out = ((eps1*(3./128))*((np.pi*Mc*x)**(-5./3))*( 1.+ (20./9)*((743./336)+(11./4)*eta)*(np.pi*Mt*x)**(2./3) -4.*(4.*np.pi - beta)*(np.pi*Mt*x) + 10.*eps2*((3058673./1016064) + (eta*5429./1008) + (617*(eta**2)/144) - sig)*(np.pi*Mt*x)**(4./3)) + 2.*np.pi*x*DT)/(2.*np.pi*factor)
    return out


# In[4]:


out_phase = f_x(xs)
plt.plot(xs,out_phase)
plt.scatter(xs,out_phase,color='k',marker='.')


# # Phase

# ## Comparing IMR Inspiral to Fergus's work

# This uses only the first 4 terms in the PN expansion, hence is not completely suitable for full IMR waveforms

# In[5]:


from numpy import pi, log
f = xsT
tc = 0.005 #0 #Time of Coalescence
p = -0.76 #-55 #Phase shift
mc = 65.0*4.925490947641266978197229498498379006e-6 #Total Mass in seconds
Mf = mc*f
phiins = 2*pi*f*tc - p + (0.192233276998395*(Mf)**0.666666666666667 - 2.19688783134229*(Mf)**1.0\
                          + 2.96024568849217*(Mf)**1.33333333333333\
                          + 0.0139119547125656)/(Mf)**1.66666666666667 - pi/4
phiins = (phiins)/(2*pi)


# In[6]:


plt.plot(xsT,phiins,color='r',label='IMR')
plt.plot(xsT,out_phase,color='k',ls='--', label='Fergus')
#plt.plot(f,out-phiins,label='Difference')
plt.legend()
plt.show()


# ## Defining Functions

# In[7]:


from numpy import pi, log, arctan, exp
fRD = 0.08803560530535628 
fDM = 0.013587622848993007
mc = 65.0*4.925490947641266978197229498498379006e-6
lim1 = 0.018/mc
lim2 = fRD/2 /mc
f = np.linspace(0.1,400.0,20000)
M = mc*f


# In[8]:


def phiins_align(m,f): #m = total mass, f = frequency array
    M = m*f
    temp = (0.192233276998395*M**0.666666666666667 - 2.19688783134229*M**1.0 + 2.96024568849217*M**1.33333333333333 + M**1.66666666666667*(14.5236327040957*log(M) + 13.7382345406982 + 14.5236327040957*log(pi)) + M**2.0*(-14.9248887589603*log(M) - 151.602171055319 - 14.9248887589603*log(pi)) + 227.418872661984*M**2.33333333333333 - 291.846189618398*M**2.66666666666667 - 2730.3817082861*M**3.0 + 3143.98200779746*M**3.33333333333333 + 9618.72578330961*M**3.66666666666667 + 0.0139119547125656)/M**1.66666666666667 - 1318.57198956936
    return temp

def phiint_align(m,f):
    M = m*f
    temp = 345.028196850348*M - 9.63599141522925*log(M) - 1416.4530276395 + 1.29058693573989e-5/M**3
    return temp

def phimrd_align(m,f):
    M = m*f
    temp = -76.9234895092469*M**0.75 + M*(68.6084707764252/fRD**0.25 + 0.94033600310695 - 19.2719828304585/fRD + 0.842931455864094/fRD**2 - 0.000619481729155146/fRD**4 - 1.06458720296323/(fDM*(1.0 + 0.243181811406409*fRD**2/fDM**2))) + 344.087860847241*M + 45.7389805176168*fRD**0.75 - fRD*(68.6084707764252/fRD**0.25 + 0.94033600310695 - 19.2719828304585/fRD + 0.842931455864094/fRD**2 - 0.000619481729155146/fRD**4 - 1.06458720296323/(fDM*(1.0 + 0.243181811406409*fRD**2/fDM**2)))/2 + 0.470168001553475*fRD - 9.63599141522925*log(fRD/2) + 1.06458720296323*arctan(0.49313467877083*fRD/fDM) + 1.06458720296323*arctan((M - 0.99313467877083*fRD)/fDM) - 1416.4530276395 - 0.421465727932047/fRD + 0.000103246954859191/fRD**3 + 0.210732863966024/M
    return temp


# In[9]:


phi_net = np.piecewise(f,[f<=lim1,f<=lim2,f>lim2],[lambda f: phiins_align(mc,f), lambda f: phiint_align(mc,f), lambda f: phimrd_align(mc,f)])
out_fin = out_phase
tc = 0.00093 #Time of Coalescence
p = -1373 #Phase shift
#These values DIFFER from the previously set values because this expression below includes all PN terms that the previous phase didn't
phi_net_plot = (2*np.pi*f*tc - np.pi/4 - p + phi_net)/(2*np.pi)


# In[10]:


plt.plot(f,phi_net_plot,label='IMR Phase',color='k')
#plt.plot(f,(phiins_align(mc,f)-p)/(2*np.pi),label='Inspiral')
#plt.plot(f,(phiint_align(mc,f)-p)/(2*np.pi),label='Intermediate')
#plt.plot(f,(phimrd_align(mc,f)-p)/(2*np.pi),label='Merger-Ringdown')
plt.plot(xsT,out_fin,label='Fergus (Only Inspiral)',color='r',linestyle='--')
plt.axvline(lim1,label='Inspiral - Intermediate',linestyle='--')
plt.axvline(lim2,label='Intermediate - Merger',linestyle='-.')
plt.legend()
plt.ylabel("$\phi (cycles)$")
plt.xlabel("f (Hz)")
plt.ylim(0,3)
plt.xlim(30,lim1)
#plt.savefig("phase2.png")
plt.show()

# In[11]:


#This is just to check whether the piecewise function has been correctly implemented
diff1 = phi_net - phiins_align(mc,f)
diff2 = phi_net - phiint_align(mc,f)
diff3 = phi_net - phimrd_align(mc,f)
plt.plot(f,diff1,label='Inspiral')
plt.plot(f,diff2,label='Intermediate',lw=3,alpha=0.5)
plt.plot(f,diff3,label='Merger-Ringdown')
plt.ylim(-2,2)
plt.legend()


# ## Choosing bins - The results from this method are not used for match calculation

# In[12]:


ph_fun = scipy.interpolate.interp1d(f,phi_net_plot) 
fmin = 40
fmax = 350
N = 500
fr = np.linspace(fmin,fmax,N)
df = (fmax-fmin)/N
der2 = scipy.misc.derivative(ph_fun,x0=fr,dx=1*df,n=2)
phases = ph_fun(fr)
plt.plot(fr,der2)
#plt.xlim(40,70)


# In[13]:


k = 0.019 #0.1 (16 bins)
bins = np.array([])
bins = np.append(bins,fmin)
i = fmin
while i<fmax:
    der2 = scipy.misc.derivative(ph_fun,x0=i,dx=3*df,n=2)
    delta = (k*abs(1/der2))**(1/2)
    fn = i + delta
    bins = np.append(bins,fn)
    i = fn


# In[14]:


print(len(bins))
bins[-1]=fmax #This is done because the final bin exceeds the maximum frequency allowed


# In[15]:


plt.plot(fr,phases)
plt.scatter(bins,ph_fun(bins),color='k',marker='.')
plt.xlabel("f (Hz)")
plt.ylabel("$\phi (cycles)$")
#plt.savefig("phase_points.png")


# In[16]:


app_ph = scipy.interpolate.InterpolatedUnivariateSpline(bins, ph_fun(bins), k=1)


# In[17]:


plt.xlabel("f (Hz)")
plt.ylabel("$\phi (cycles)$")
plt.plot(fr, app_ph(fr) - phases)
plt.title('Difference between the spline and the IMR phase')
#plt.savefig("phase_errors.png")


# ## Alternate Method - Used in match calculation

# In[18]:


ph_fun_alt = scipy.interpolate.InterpolatedUnivariateSpline(f,phi_net_plot,k=5) 
fmin = 40
fmax = 350
N = 500
fr = np.linspace(fmin,fmax,N)
df = (fmax-fmin)/N
der2_alt = ph_fun_alt.derivative(n=2)
phases = ph_fun_alt(fr)
der = der2_alt(fr)
#plt.plot(fr,der)
#plt.plot(fr,der2)
#plt.xlim(40,70)


# In[19]:


k = 0.018 #0.1 (16 bins)
bins = np.array([])
bins = np.append(bins,fmin)
i = fmin
while i<fmax:
    der2 = der2_alt(i)
    delta = (k*abs(1/der2))**(1/2)
    fn = i + delta
    bins = np.append(bins,fn)
    i = fn


# We need 33 points for 32 bins

# In[20]:


print(len(bins))
bins[-1] = fmax


# In[21]:


app_ph = scipy.interpolate.InterpolatedUnivariateSpline(bins, ph_fun_alt(bins), k=1)


# In[22]:


plt.xlabel("f (Hz)")
plt.ylabel("$\phi (cycles)$")
plt.plot(fr, app_ph(fr) - phases)
plt.title('Difference between the spline and the IMR phase')


# ### Save data

# In[23]:


#RUN THIS TO SAVE DATA LOCALLY

#np.savetxt('Phase_bins.txt',bins)
#np.savetxt('Phase_bins_phases.txt',ph_fun(bins))

print("AMPLITUDE")
# # Amplitude

# ## Defining Functions

# In[24]:


def ampins(m,f):
    M = m*f
    temp = 1.0*(-1.65346016027308*M**0.666666666666667 - 20.9165547229562*M**1.33333333333333 - 49.7411328692393*M**2.0 + 3732.32311535403*M**2.33333333333333 - 38985.411777689*M**2.66666666666667 + 97150.3925991756*M**3.0 + 1)/M**1.16666666666667
    return temp

def ampint(m,f):
    M = m*f
    temp = 1.0*(-24857.8876111804*M**4 - 555.774331218285*M**3 + 527.285126438352*M**2 - 24.8870827483312*M + 1.07021302528877)/M**1.16666666666667
    return temp

def ampmrd(m,f):
    M = m*f
    temp = 0.0223314133093818*fDM*exp(-0.644487671942311*(M - fRD)/fDM)/(M**1.16666666666667*(1.74724436683344*fDM**2 + (M - fRD)**2))
    return temp


# In[25]:


lim1a =  0.014/mc #Inspiral - Intermediate boundary
lim2a =  0.07803912805383108/mc #Intermediate - Merger boundary
amp0 = 0.33733861901045753
plt.show()

# In[26]:


amp_net = amp0*np.piecewise(f,[f<=lim1a,f<=lim2a,f>lim2a],[lambda f: ampins(mc,f), lambda f: ampint(mc,f), lambda f: ampmrd(mc,f)])


# In[27]:


plt.loglog(f,amp_net/(45/0.25),label='IMR Amplitude',color='k')
plt.loglog(f,amp_net,label='IMR Amplitude',color='red')
plt.axvline(lim1a,label='Inspiral - Intermediate',linestyle='--')
plt.axvline(lim2a,label='Intermediate - Merger',linestyle='-.')
plt.legend()
plt.xlabel('f (Hz)')
plt.ylabel('A')
plt.show()
#plt.xlim(40,200)
#plt.savefig("amp2.png")
fr = np.linspace(40,350,311)
amp_net2 = amp0*np.piecewise(fr,[fr<=lim1a,fr<=lim2a,fr>lim2a],[lambda fr: ampins(mc,fr), lambda fr: ampint(mc,fr), lambda fr: ampmrd(mc,fr)])
plt.loglog(fr,amp_net2,label='IMR Amplitude',color='red')
plt.show()
print('frequency')
print(fr)
results = []
for i,value in enumerate(fr):
    temp = [value,amp_net2[i]]
    results.append(temp)
with open("test", "wb") as fp:
    pickle.dump(results, fp)
# ## Choosing bins - The results from this method are not used for match calculation

# In[28]:


#The second derivative of the amplitude has a discontinuity at the Intermediate - Merger boundary.
am_fun = scipy.interpolate.interp1d(f,amp_net)
fmin = 40
fmax = 350
N = 500
fr = np.linspace(fmin,fmax,N)
df = (fmax-fmin)/N
der2_am = scipy.misc.derivative(am_fun,x0=fr,dx=1*df,n=2,)
amps = am_fun(fr)
plt.plot(fr,der2_am,color='k')
plt.axvline(lim2a,label='Intermediate - Merger',linestyle='-.')
#plt.plot(fr,der2)
#plt.xlim(40,70)


# In[29]:


k = 0.2708
bins_amp = np.array([])
bins_amp = np.append(bins_amp,fmin)
i = fmin
while i<fmax:
    der2 = scipy.misc.derivative(am_fun,x0=i,dx=3*df,n=2)
    delta = (k*abs(1/der2))**(1/2)
    fn = i + delta
    bins_amp = np.append(bins_amp,fn)
    i = fn


# In[30]:


# I had to manually set some of the bins due to the discontinuity in the derivative, which isn't ideal. I am trying to fix this
print(len(bins_amp))
bins_amp[24] = 210.0
bins_amp[28] = 300.0
bins_amp[29] = 319.95
bins_amp[30] = 332.495
bins_amp[31] = 347.61
print(bins_amp)


# In[31]:


plt.loglog(fr,amps)
plt.scatter(bins_amp,am_fun(bins_amp),color='k',marker='.')
plt.xlabel("f (Hz)")
plt.ylabel("$A$")
#plt.savefig("Amplitude_points.png")


# In[32]:


app_amp = scipy.interpolate.InterpolatedUnivariateSpline(bins_amp, am_fun(bins_amp), k=1)


# In[33]:


plt.xlabel("f (Hz)")
plt.ylabel("$A$")
plt.plot(fr, app_amp(fr) - amps)
#plt.savefig("amp_errors.png")


# ## Alternate method - Used in match calculation

# In[34]:


am_fun_alt = scipy.interpolate.UnivariateSpline(f,amp_net,k=5,s=0.003) 
fmin = 40
fmax = 350
N = 500
fr = np.linspace(fmin,fmax,N)
df = (fmax-fmin)/N
der2_alt = am_fun_alt.derivative(n=2)
amps = am_fun_alt(fr)
der = der2_alt(fr)
#plt.loglog(fr,amps)
#plt.loglog(f,amp_net)
#plt.xlim(fmin,fmax)
#plt.ylim(am_fun_alt(fmax)-0.1,am_fun_alt(fmin)+1)
plt.axvline(lim2a,color='k',alpha=0.5)
plt.plot(fr,der)
n = 0.1
#plt.ylim(-n/4,n)


# In[35]:


#k = 0.2705
k=0.001
bins_amp = np.array([])
bins_amp = np.append(bins_amp,fmin)
i = fmin
while i<fmax:
    der2 = der2_alt(i)
    delta = (k*abs(1/der2))**(1/2)
    fn = i + delta
    bins_amp = np.append(bins_amp,fn)
    i = fn


# Again, we need 33 points for 32 bins

# In[36]:


print(len(bins_amp))
bins_amp[-1] = fmax
bins_amp


# In[37]:


fig,ax = plt.subplots(1)
ax.loglog(fr,amps,color='k',label='IMR Amplitude')
ax.axvline(lim1a,alpha=0.5,ls='--')
ax.axvline(lim2a,alpha=0.5,ls='--')
ax.fill([fmin-5,lim1a,lim1a,fmin-5],[0,0,47.0,47.0],color='gray',alpha=0.5,label='Inspiral')
ax.fill([lim1a,lim2a,lim2a,lim1a],[0,0,47.0,47.0],color='orange',alpha=0.5,label='Intermediate')
ax.fill([lim2a,fmax+5,fmax+5,lim2a],[0,0,47.0,47.0],color='red',alpha=0.5,label='Merger-Ringdown')
plt.scatter(bins_amp,am_fun_alt(bins_amp),color='k',marker='.')
plt.xlim(fmin-5,fmax+5)
plt.ylim(0,47.0)
plt.legend()
plt.xlabel("f (Hz)")
plt.ylabel("$A$")
#plt.savefig('IMRAmp.png')
plt.show()
data=[]
print(am_fun_alt(bins_amp))
amplitude=am_fun_alt(bins_amp)
for i,value in enumerate(bins_amp):
    data.append([value,amplitude[i]])
print(data)
# In[38]:


app_amp = scipy.interpolate.InterpolatedUnivariateSpline(bins_amp, am_fun_alt(bins_amp), k=1)


# In[39]:


plt.xlabel("f (Hz)")
plt.ylabel("$A$")
plt.plot(fr, app_amp(fr) - amps)


# ### Save data

# In[40]:


#RUN THIS TO SAVE DATA LOCALLY

#np.savetxt('Amp_bins.txt',bins_amp)
#np.savetxt('Amp_bins_amplitudes.txt',am_fun(bins_amp))


# # Building and checking waveforms

# In[41]:


ref_freq = 40


# This is used to set the waveforms to zero below a particular reference frequency. 
# This is necessary because the equations we use to calculate the match (as below) involve a fourier transform.
# 
# $\mathcal{O}(h_1,h_2) = \frac{4}{||h_1|| ||h_2||} \max_{t_0}\left|\mathcal{F}^{-1}\left[\frac{\tilde{h}_1^*(f) \tilde{h}_2(f)}{S_n(f)}\right](t_0)\right|$

# In[42]:


def model(fre): #IMRPhenomD model for the waveform
    if fre < ref_freq:
        return 0
    temp = am_fun_alt(fre)*np.exp(1j*m*ph_fun_alt(fre))
    return temp
def approx(fre): #Linear Spline approximant for the waveform
    if fre < ref_freq:
        return 0
    temp = app_amp(fre)*np.exp(1j*m*ph_fun_alt(fre))
    return temp


# In[43]:


f = np.linspace(0,350,2000)
m = 1
waveform_model = np.array([])
waveform_approx = np.array([])
for fi in f:
    tm = model(fi)#am_fun_alt(i)*complex(np.cos(m*ph_fun_alt(i)),np.sin(m*ph_fun_alt(i)))
    ta = approx(fi)#app_amp(i)*complex(np.cos(m*app_ph(i)),np.sin(m*app_ph(i)))
    waveform_model = np.append(waveform_model,tm)
    waveform_approx = np.append(waveform_approx,ta)


# In[44]:


mod_plot = np.real(waveform_model)
app_plot = np.real(waveform_approx)
mod_iplot = np.imag(waveform_model)
app_iplot = np.imag(waveform_approx)
#plt.plot(f,mod_plot,color='k',label='Model')
#plt.plot(f,app_plot,color='r',label='Approximant',ls='--')
plt.plot(f,mod_plot,color='k',label='Model')
plt.plot(f,app_plot,color='r',label='Approximant',ls='--')
#plt.plot(f, mod_plot-app_plot)
#plt.plot(f, mod_iplot-app_iplot)
plt.xlim(0,350)
#plt.ylim(-5,27)
plt.legend()


# # Mismatch Calculation

# ### PSD

# In[45]:


#psd = np.genfromtxt("GWTC1_GW150914_PSDs.dat") #LIGO PSD
psd = np.genfromtxt("aLIGODesign.txt") #Advanced LIGO PSD - This is what we are using
psd = np.transpose(psd)


# In[46]:


freqs = psd[0]
psd_arr = psd[1]


# In[47]:


psdfun = scipy.interpolate.interp1d(freqs,psd_arr,fill_value='extrapolate') 
f = np.linspace(0, 350, 10000)
plt.loglog(f,psdfun(f))


# ### Match and Mismatch

# In[48]:


h1 = np.array([approx(fi) for fi in f])
h2 = np.array([model(fi) for fi in f])
psd_arr = psdfun(f)


# We also calculated the matches with a padded integrand to check if we got more accurate results. The results don't vary by much in the end.

# In[49]:


integrand = h1*h2.conj()/psd_arr
h1_norm = np.abs(4*np.dot(h1.conj(), h1/psd_arr))
h2_norm = np.abs(4*np.dot(h2.conj(), h2/psd_arr))

p = 4*np.abs(np.fft.fft(integrand)) / np.sqrt(h1_norm*h2_norm)
integrand_padded = np.lib.pad(integrand, 2*len(integrand))
p_zpad = 4*np.abs(np.fft.fft(integrand_padded)) / np.sqrt(h1_norm*h2_norm)


# In[50]:


h1_norm, h2_norm


# In[51]:


match, matchpad = np.max(p), np.max(p_zpad)
print(match, matchpad)


# In[52]:


mismatch = (1 - match)*100
mismatchpad = (1 - matchpad)*100
print(mismatch,mismatchpad) #In Percentages


# In[53]:


print(np.argmax(p_zpad))
plt.plot(p);
plt.plot(p_zpad);


# This implies that the time shift that gives us the maximum match is 0.

# # Calculating $\theta$ for the amplitudes
# This makes use of Eqn (5). from Fergus' paper with a small modification. The probabilities are calculated  by integrating the linear spline approximant between the chosen frequencies. Alternatively, we can calculate the probabilities using the area of a trapezium. These probabilities are then added according to the limits in the equation. The notation in the following block has been chosen to match the notation used in the paper where possible. 

# In[54]:


n = 5
m = np.array([0,1,2,3,4])
costhetamj = []

for mi in m:
    j = np.arange(0,((2**mi)),1)
    temp = np.array([])
    for ji in j:
        
        low = int(ji*(2**(n-mi)))
        high_num = int((ji+0.5)*(2**(n-mi)))
        high_den = int((ji+1)*(2**(n-mi)))
        
        area_num = [(1/2)*(bins_amp[i+1] - bins_amp[i])*(app_amp(bins_amp[i+1]) + app_amp(bins_amp[i])) for i in np.arange(low,high_num,1)]
        area_den = [(1/2)*(bins_amp[l+1] - bins_amp[l])*(app_amp(bins_amp[l+1]) + app_amp(bins_amp[l])) for l in np.arange(low,high_den,1)]
        
        cos2th = np.sum(area_num)/np.sum(area_den)
        costh = np.sqrt(cos2th)
        th = np.arccos(costh)
        
        dict = {"m":mi,"j":ji,"theta":th}
        temp = np.append(temp,dict)
    
    costhetamj = np.append(costhetamj,temp)


# In[55]:


costhetamj

