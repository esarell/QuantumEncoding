from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, Aer, IBMQ
from scipy.interpolate import approximate_taylor_polynomial
from qiskit.circuit.library import RGQFTMultiplier, DraperQFTAdder, ExactReciprocal
from qiskit.circuit.library.basis_change import QFT as QFT_pre
from qiskit.extensions import HamiltonianGate
from qiskit.circuit.library.standard_gates import PhaseGate, RYGate, CSwapGate
from qiskit.circuit.library.arithmetic.integer_comparator import IntegerComparator
from qiskit.circuit.library.standard_gates import XGate
from qiskit.quantum_info import Pauli
from qiskit.circuit.library.boolean_logic import OR
from qiskit.quantum_info import random_unitary
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Unroller
import itertools as it
import numpy as np
import sympy
import itertools as it
import matplotlib.pyplot as plt
#from LPF_coefficients import get_bound_coeffs

'''
The Grover-Rudolph code works by first calling the sim_Af_GR() function which sets up the distribution,
we will then encode. Then it calls the Grover_Rudolph_func() function.
This then calls the piecewise_function_posmulti() function
Which then calls the label_gate() function. Which is made up of the integer_compare() and increment_gate()
function. 

This file also contains a few other functions to help:
QFT
bin_to_dec
QFTBinaryAdd
my_binary_repr
twos_compliment
'''

def QFT(circ, qreg, do_swaps=True, approximation_degree=0., insert_barriers=False, wrap=False, inverse=False, label='QFT'):

    n = len(qreg)

    if inverse:
        wrap = True

    if wrap:
        qreg = QuantumRegister(n, 'q_reg')
        circ = QuantumCircuit(qreg)

    circ.append(QFT_pre(n, do_swaps=do_swaps, approximation_degree=approximation_degree, insert_barriers=insert_barriers).to_gate(), qreg)

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'

    return circ

def bin_to_dec(binary, nint=None, phase=False, signmag=False):
    """
    Convert a binary string to a floating point number
    Using twos complement
    binary - input binary (string)
    nint - number of integer bits. Default to all (int)
    phase = is if its a negative balue
    """

    basis = 0. 

    n = len(binary)
    
    if phase:
        #if binary[0]=='1':
        #    sign = -1.
        #elif binary[0]=='0':
        #    sign = 1.
        if binary[0]=='1':
            if nint is None:
                nint_ = n-1
            else:
                nint_ = nint
            basis = -(2.**(nint_))
        binary = binary[1:]
        #if nint is not  None:
        #    nint += 1

    n = len(binary)
    if nint is None:
        nint = n

    digit = 0.
    for i,bit in enumerate(np.arange(nint-n,nint)[::-1]):#enumerate(np.arange(nint-n+1,nint+1)[::-1]):
        digit+=(2.**bit)*int(binary[i])
        #print(i,bit,2.**bit,binary[i],digit)
    return digit + basis

def QFTBinaryAdd(circ, qreg, binary, wrap=False, inverse=False, QFT_on=True, iQFT_on=True, label='AddA'):
    """
    |qreg> -> |qreg + A>
    """
    n1 = len(qreg)
    n2 = len(binary)

    if inverse:
        wrap = True

    if wrap:
        qreg = QuantumRegister(n1, 'q_reg1')
        circ = QuantumCircuit(qreg)

    if QFT_on:
        circ.append(QFT(circ, qreg, do_swaps=False, wrap=True), qreg[:])

    for i in np.arange(1, n1 + 1):
        for k in np.arange(1, n2 + 1):
            lam = (2. * np.pi) / (2. ** (i + k - n2))
            if lam%(2*np.pi)==0.:
                continue
            if binary[i-1]=='1':
                circ.p(lam,qreg[k-1])

    if iQFT_on:
        circ.append(QFT(circ, qreg, do_swaps=False, wrap=True, inverse=True), qreg[:])

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'

    return circ

def my_binary_repr(digit, n, nint=None, phase=False, nround=True, overflow_error=True):
    """
    Convert a floating point digit to binary string
    in the format of twos complement
    digit - input number (float)
    n - total number of bits (int)
    nint - number of integer bits. Default to lowest required (int)
    """

    if nint is None:# or nint==n:
        if phase:
            nint = n - 1
        else:
            nint = n

    if phase:
        p = n - nint - 1
        dmax = 2.**(nint) - 2.**(-p)
        dmin = -2.**(nint)
    else:
        p = n - nint
        dmax = 2.**(nint) - 2.**(-p)
        dmin = 0.

    if overflow_error:
        if digit>dmax or digit<dmin:
            raise ValueError('Digit '+str(digit)+' does not lie in the range:',dmin,'-',dmax,n,nint,p)

    if nround:
        n += 1
        p += 1

    value = digit
    bin_out = ''
    if phase:
        if value<0.:
            value+=2.**nint
            bin_out+='1'
        else:
            bin_out+='0'
    
    for i,bit in enumerate(np.arange(-p,nint)[::-1]):
        bin_out+=str(int(np.floor(value/2.**bit)))
        if value>=2.**bit:
            value-=2.**bit

    if nround:
        carry = True
        bin_out = np.array(list(bin_out))
        for i in np.arange(n)[::-1]:
            if not carry:
                break
            if bin_out[i]=='1':
                bin_out[i]='0'
            elif bin_out[i]=='0':
                bin_out[i]='1'
                carry = False
        bin_out = ("").join(list(bin_out[:-1]))

    return bin_out

def twos_compliment(binary):
    n = len(binary)

    if np.sum(np.array(list(binary)).astype(int))==0:
        compliment = binary
    else:
        compliment = my_binary_repr(bin_to_dec(''.join(list(np.logical_not(np.array(list(binary)).astype(bool)).astype(int).astype(str))), nint=None, phase=False) + 1, n, nint=None, phase=False, overflow_error=False)
    print(compliment)
    return compliment

def increment_gate(circ, qreg, wrap=False, inverse=False, QFT_on=True, iQFT_on=True, ncut=0, label='inc'):
    #qreg here actually refers to qlab our label gate 
    #We are currently in the phase space cause we have already applied the QFT
    n = len(qreg)

    if inverse:
        wrap = True

    if wrap:
        qreg = QuantumRegister(n, 'q_reg')
        circ = QuantumCircuit(qreg)
    
    if n<ncut:
        for i in np.arange(n)[1:][::-1]:
            if i!=0:
                xgate = XGate().control(i)
                circ.append(xgate, [*qreg[:i+1]])
    
        circ.x(qreg[0])

    else:
        '''Currently this bit of code runs every time'''
        #Adds one to the current value in the qlab
        bin_one = ''.join(np.zeros(n-1).astype(int).astype(str))+'1'
        inc_gate = QFTBinaryAdd(circ, qreg, bin_one, wrap=True, QFT_on=QFT_on, iQFT_on=iQFT_on)
        circ.append(inc_gate, qreg)

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'

    return circ

def integer_compare(circ, qreg, qtarg, qans, value, wrap=True, inverse=False, uncomp=True, label='P'):
    
    n = len(qreg)

    if wrap:
        qreg = QuantumRegister(n, 'q_reg')
        qans = QuantumRegister(n, 'q_ans')
        qtarg = QuantumRegister(1,'qtarg')
        circ = QuantumCircuit(qreg, qans,qtarg)
    
    if value<=0.:
        #geq is currently set to False
        if True:
            circ.x(qtarg)
    elif value < np.power(2,n):
        if n>1:
            #Does twos compliments of the bound we have sent in
            #returns it backwards so the first bit represents the lowest value
            twos = np.array(list(twos_compliment(my_binary_repr(value, n=n, phase=False))))[::-1]
            #goes through the binary size of the list e.g. 0,1,2,3,4,5
            for i in np.arange(n):
                if i==0:
                    #Checks if the first bit is a one
                    if twos[i]=='1':
                        #puts a 1 in the ancillary register depended on the value in the qreg
                        #ancillary keeps track of if two 1's are present and we end up with a carry
                        circ.cx(qreg[i], qans[i])
                elif i<n-1:
                    if twos[i]=='1':
                        #Check if there is a carry already present in the previous ancillary register
                        #Or if there if the value of x at this point is also a one
                        circ.compose(OR(2), [qreg[i], qans[i-1], qans[i]], inplace=True);
                    else:
                        #If twos compliment at this point is 0 then we check if both the x at this point 
                        #and previson ancillary is a 1 if so then we put a 1 in the current ancillary regsiter
                        circ.ccx(qreg[i], qans[i - 1], qans[i]);
                #Once we have reached the last bit we do the same checks but now but the carry bit
                #In the target register instead of the ancilary register
                else:
                    if twos[i]=='1':
                        circ.compose(OR(2), [qreg[i], qans[i-1], qtarg[0]], inplace=True);
                    else:
                        circ.ccx(qreg[i], qans[i - 1], qtarg[0]);
            #if not geq:
            #    circ.x(qtarg);
            
            if uncomp:
                for i in np.arange(n-1)[::-1]:
                    if i==0:
                        if twos[i]=='1':
                            circ.cx(qreg[i], qans[i]);
                    else:
                        if twos[i]=='1':
                            circ.compose(OR(2), [qreg[i], qans[i-1], qans[i]], inplace=True);
                        else:
                            circ.ccx(qreg[i], qans[i-1], qans[i]);
        else:
            circ.cx(qreg[0], qtarg)
        
            #if not geq:
                #circ.x(qtarg);
    #else:
        #if not geq:
            #circ.x(qtarg);

    '''if wrap:
        circ = circ.to_gate()
        circ.label = label'''

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'

    return circ


def label_gate(circ, qreg, qtarg, qans, qlab, bounds=None, wrap=False, nint=None, inverse=False, phase=False, ncut=0, label='LABEL'):
    n = len(qreg)

    if nint is None:
        nint = n
        if phase:
            nint -= 1

    if inverse:
        wrap = True

    if bounds is None:
        bounds = [0, 2**n]
    print("bounds_label:",bounds)
    nlab = int(np.ceil(np.log2(len(bounds))))
    print("nlab:",nlab)
    if nlab!=qlab.size:
        raise ValueError('Size of label register does not match the number of bounds placed.')

    if len(qans)!= n+1:
        raise ValueError('Ancilla register must have one more qubit than input register.')

    if wrap:
        qreg = QuantumRegister(n, 'q_reg')
        qtarg = QuantumRegister(1, 'q_targ')
        qans = QuantumRegister(n+1, 'q_ans')
        qlab = QuantumRegister(nlab, 'q_lab')
        circ = QuantumCircuit(qreg, qtarg, qans, qlab)
  
    circ.x(qreg[-1]);

    qreg = [*qreg, qans[-1]]
    qans = [*qans[:-1]]

    
    if nlab>=ncut:
        circ.append(QFT(circ, qlab, do_swaps=False, wrap=True), qlab[:])

    for i,bound_ in enumerate(bounds):
        #Flips the first bit so 1->9 ect up to 8 then 8->0 up to 15 for m =4
        #Current working theory is that this means we only label the second half of
        #the distibution
        print("indivdual bound:",bound_)
        binary = my_binary_repr(bound_, n=n, nint=nint, phase=phase)
        print("binary before",binary)
        if binary[0]=='0':
            binary = '1'+binary[1:]
            print("Binary bit chnage:",binary)
        elif binary[0]=='1':
            binary = '0'+binary[1:]
            print("Binary bit chnage:",binary)

        #Converts back to binary
        bound = bin_to_dec(binary, nint=None, phase=False)
        print("ind_bound_change:",bound)
        intcomp_gate = integer_compare(circ, qreg, qtarg, qans, bound, geq=True, wrap=True, uncomp=False, label='P'+str(i))
        circ.append(intcomp_gate, [*qreg, qtarg[0], *qans[:]])

        inc_gate = increment_gate(circ, qlab, wrap=True, label='SET'+str(i), ncut=ncut, QFT_on=False, iQFT_on=False).control(1)
        circ.append(inc_gate, [qtarg[0], *qlab[:]])

        intcomp_gate_inv = integer_compare(circ, qreg, qtarg, qans, bound, geq=True, wrap=True, uncomp=False, inverse=True, label='P'+str(i))
        circ.append(intcomp_gate_inv, [*qreg, qtarg[0], *qans[:]])

    if nlab>=ncut:
        circ.append(QFT(circ, qlab, do_swaps=False, wrap=True, inverse=True), qlab[:])

    circ.x(qreg[-2]);

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'

    
    return circ

def piecewise_function_posmulti(circ, q_x, q_y, q_lab, q_coff, coeffs, bounds, nint=None, nintx=None, nintcs=None, phase=False, wrap=False, inverse=False, unlabel=False, unfirst=False, comp2=False, label='f_x'):

    nlab = int(np.ceil(np.log2(len(bounds))))
    nx = len(q_x)
    n = len(q_y)
    nc = len(q_coff)

    Nord, Ncoeffs = coeffs.shape

    if Nord!=2:
        raise ValueError('Currently only working for linear piecewise approximations.')

    if nint is None:
        nint = n

    if nlab!=len(q_lab):
        raise ValueError('Size of label register is smaller than number of bounds.')
        return 0

    if np.any(nintcs==None):
        nintcs = []
        for coeffs_ in coeffs:
            nintcs.append(int(np.ceil(np.log2(np.max(np.abs(coeffs_))))))
        nintcs[-1] = nint
        nintcs = np.array(nintcs).astype(int)

    if nintx is None:
        nintx = int(np.ceil(np.log2(np.max(np.abs(bounds)))))

    if inverse:
        wrap = True

    if wrap:
        q_x = QuantumRegister(nx, 'q_x')
        q_y = QuantumRegister(n, 'q_y0')
        q_lab = QuantumRegister(nlab, 'q_lab')
        q_coff = QuantumRegister(nc, 'q_coff')

        circ = QuantumCircuit(q_x, q_y, q_lab, q_coff)

    q_ans = [*q_coff, *q_y]

    if len(q_ans)<nx:
        raise ValueError('Coefficient/output register must be greater than half the length of the x register.')
    print("Bounds before label:",bounds)
    lab_gate = label_gate(circ, q_x, q_ans[0], q_ans[1:nx+2], q_lab, bounds=bounds, nint=nintx, phase=False, wrap=True)
    circ.append(lab_gate, [*q_x, *q_ans[:nx+2], *q_lab]);

    #Currently commented all this out as im just trying to figrue out the label gate
    '''
    #y1in_gate = first_gate(circ, q_y, q_lab, coeffs[1], nint=nintcs[0,1], phase=phase, wrap=True)
    #circ.append(y1in_gate, [*q_y, *q_lab]);

    #circ.append(QFT(circ, q_y, do_swaps=False, wrap=True), q_y);
    
    for i in np.arange(1, Nord):

        y0in_gate = first_gate(circ, q_coff, q_lab, coeffs[i-1], nint=nintcs[i-1,0], phase=phase, wrap=True, comp2=comp2)
        circ.append(y0in_gate, [*q_coff, *q_lab]);

        mul_gate = QFTPosMultiplicand(circ, q_coff, q_x, q_y, wrap=True, nint1=nintcs[i-1,0], nint2=nintx, nint3=nint, iQFT_on=True, QFT_on=True, comp2=comp2)
        circ.append(mul_gate, [*q_coff, *q_x, *q_y]);

        if unfirst or True:
            y0in_gate_inv = first_gate(circ, q_coff, q_lab, coeffs[i-1], nint=nintcs[i-1,0], phase=phase, wrap=True, comp2=comp2, inverse=True)
            circ.append(y0in_gate_inv, [*q_coff, *q_lab]);

        y1in_gate = first_gate(circ, q_coff, q_lab, coeffs[i], nint=nintcs[i-1,1], phase=phase, wrap=True)
        circ.append(y1in_gate, [*q_coff, *q_lab]);

        add_gate = QFTAddition(circ, q_coff, q_y, wrap=True, phase=phase, nint1=nintcs[i-1,1], nint2=nint, QFT_on=True, iQFT_on=True)
        circ.append(add_gate, [*q_coff, *q_y]);
    '''
    #circ.append(QFT(circ, q_y, do_swaps=False, wrap=True, inverse=True), q_y);

    if unlabel:
        lab_gate_inv = label_gate(circ, q_x, q_ans[0], q_ans[1:nx+2], q_lab, bounds=bounds, nint=nintx, phase=False, wrap=True, inverse=True)
        circ.append(lab_gate_inv, [*q_x, *q_ans[:nx+2], *q_lab]);

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'

    return circ




def Grover_Rudolph_func(circ, qx, qanc, qlab, qcoff, probs, wrap=False, inverse=False, mtol=1, mmax=None, norder=1, label='GR_func'):

    nx = len(qx)
    n = len(qcoff)
    nanc = len(qanc)
    nlab = len(qlab)
    print("nx:",nx)
    print("Nlan!",nlab)
    if mmax==None:
        mmax=nx

    if n!=nanc:
        raise ValueError('I think ancilla and coefficient reg should be the same size')
    
    if len(probs)!=2**nx:
        raise ValueError('Probabilities must equal length of x register')

    if inverse:
        wrap = True

    if wrap:
        qx = QuantumRegister(nx, 'q_x')
        qanc = QuantumRegister(nanc, 'q_anc')
        qlab = QuantumRegister(nlab, 'q_lab')
        qcoff = QuantumRegister(nanc, 'q_coff')
        circ = QuantumCircuit(qx, qanc, qlab, qcoff)


    print("mtol,",mtol)
    print("mmax",mmax)
    for m in np.arange(mtol,mmax):
        def GR_func(j):
            j = np.array(j).astype(int)
            As = []
            for i in np.arange(2**m):
                As.append(np.sum(probs[i*2**(nx-m):(i+1)*2**(nx-m)]))
            As1 = []
            for i in np.arange(2**(m+1)):
                As1.append(np.sum(probs[i*2**(nx-(m+1)):(i+1)*2**(nx-(m+1))]))
            return np.arccos(np.sqrt(np.array(As1)[::2][j]/np.array(As)[j]))

        js = np.arange(2**m)
        coeffs = GR_func(js)
        print("m:",m)
        #bounds_ = np.linspace(0,2**m,(2**nlab)+1).astype(int)
        bounds_ = np.linspace(0,2**m-1,(2**nlab)).astype(int)
        coeffs = get_bound_coeffs(GR_func, bounds_, norder, reterr=False).T#[::-1]
        bounds = bounds_[1:]
        print("bounds_GR:",bounds)
        max_list0 = np.array([coeffs[0], coeffs[1], coeffs[0]*2**nx, (coeffs[0]*2**nx)+coeffs[-1]])
        max_list1 = max_list0
        nintcs = []
        nintcs.append(int(np.ceil(np.log2(np.max(np.abs(max_list0))))))
        nintcs.append(int(np.ceil(np.log2(np.max(np.abs(max_list1))))))
        nint = nintcs[-1]
        nintcs = np.array([nintcs])
        
        func_gate = piecewise_function_posmulti(circ, qx[nx-m:], qanc, qlab, qcoff, coeffs, bounds, nint=nint, nintx=nx, nintcs=nintcs, phase=False, wrap=True, unlabel=False, unfirst=False)
        circ.append(func_gate, [*qx[nx-m:], *qanc, *qlab, *qcoff])

        func_gate_inv = piecewise_function_posmulti(circ, qx[nx-m:], qanc, qlab, qcoff, coeffs, bounds, nint=nint, nintx=nx, nintcs=nintcs, phase=False, wrap=True, inverse=True, unlabel=False, unfirst=False)
        circ.append(func_gate_inv, [*qx[nx-m:], *qanc, *qlab, *qcoff])

    for m in np.arange(mmax,nx):
        circ.h(qx[nx-m-1]);

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'
    circ.draw('mpl')
    plt.plot()
    return circ

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
width=0.75
color='black'
fontsize=28
ticksize=22
figsize=(10,10)


def sim_Af_GR(nx, nlab, ncut0, ncut1, nsig0, nsig1, statename='./state_vectors/psi_f', fmin=40., fmax=168., m1=35., m2=30., Tfrac=100., beta=0., sig=0.):

    phase = True

    nintx = int(np.ceil(np.log2((fmax-fmin))))
    print("nintx:",nintx)

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

    def amplitude(nqubit):
        xs = np.linspace(fmin, fmax, 2**nqubit)
        amps = xs**(-7./6)
        norm = np.sqrt(np.sum(np.abs(amps)**2))
        return amps/norm

    def f_x(x, eps1=1, eps2=1, factor=1):
        x = x_trans(x)
        out = ((eps1*(3./128))*((np.pi*Mc*x)**(-5./3))*( 1.+ (20./9)*((743./336)+(11./4)*eta)*(np.pi*Mt*x)**(2./3) -4.*(4.*np.pi - beta)*(np.pi*Mt*x) + 10.*eps2*((3058673./1016064) + (eta*5429./1008) + (617*(eta**2)/144) - sig)*(np.pi*Mt*x)**(4./3)) + 2.*np.pi*x*DT)/(2.*np.pi*factor)
        return out
    print('Qubits before:', nx, nlab)
    n, nc, nlab, nint_unused, nintcs_unused, coeffs_unused, bounds_unused = optimize_coeffs_qubits(f_x, nx, nlab, nintx, ncut0, ncut1, nsig0=nsig0, nsig1=nsig1)
    #print("bounds",bounds)
    print('Qubits:', nx, n, nc, nlab)
    #print('Integer qubits:', nintx, nint, nintcs[0,0], nintcs[0,1])
    print('Memory:', 16*(2**(nc+n+nx+nlab))/2**20)

    if 16*(2**(nc+n+nx+nlab))/2**20>7568:
        raise ValueError('Too many qubits!',nc+n+nx+nlab)

    #####################################################################

    q_x = QuantumRegister(nx, 'q_x')
    q_y = QuantumRegister(n, 'q_y')
    q_lab = QuantumRegister(nlab, 'q_lab')
    q_coff = QuantumRegister(nc, 'q_coff')
    q_classical = ClassicalRegister(1,"classical")
    circ = QuantumCircuit(q_x, q_y, q_lab, q_coff,q_classical)

    q_ycoff = [*q_y, *q_coff]

    GR_gate = Grover_Rudolph_func(circ, q_x, q_ycoff[:(nc+n)//2], q_lab, q_ycoff[(nc+n)//2:2*((nc+n)//2)], amplitude(nx)**2,q_classical ,mtol=4, mmax=6, wrap=True)
    circ.append(GR_gate, [*q_x, *q_ycoff[:(nc+n)//2], *q_lab, *q_ycoff[(nc+n)//2:2*((nc+n)//2)]]);

    if list(dict(circ.decompose(reps=10).count_ops()).keys())!=['u', 'cx']:
        raise ValueError('Cannot decompose circuit into u and CX gates.')

    print('CNOT gate count:',dict(circ.decompose(reps=10).count_ops())['cx'])

    backend = Aer.get_backend('statevector_simulator')
    job = execute(circ, backend)
    result = job.result()
    state_vector = result.get_statevector()

    state_vector = np.asarray(state_vector).reshape((2**nc,2**nlab,2**n,2**nx))
    state_v = state_vector[0,0,0,:].flatten()
    print(np.sum(np.abs(state_v)**2))

    np.save(statename, state_v)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Plot the simulated state vector compared to target state")
    parser.add_argument('--statename', help="Numpy file with saved state vector.", default='./state_vectors/amp_state')
    parser.add_argument('--nx', help="Number of frequency qubits.", default=6)
    parser.add_argument('--nlab', help="Number of label qubits.", default=4)
    parser.add_argument('--ncut0', help="Number of bits to round coefficient 0 to.", default=7)
    parser.add_argument('--ncut1', help="Number of bits to round coefficient 1 to.", default=64)
    parser.add_argument('--nsig0', help="Round coefficient 0 to given significant figure.", default=10)
    parser.add_argument('--nsig1', help="Round coefficient 1 to given significant figure.", default=10)
    parser.add_argument('--fmin', help="Minimum frequency.", default=40.)
    parser.add_argument('--fmax', help="Maximum frequency.", default=168.)
    parser.add_argument('--m1', help="Component mass 1.", default=30.)
    parser.add_argument('--m2', help="Component mass 2.", default=35.)
    parser.add_argument('--beta', help="Beta spin parameter.", default=0.)
    parser.add_argument('--sig', help="Sigma spin parameter.", default=0.)
    parser.add_argument('--Tfrac', help=" ", default=100.)

    opt = parser.parse_args()

    sim_Af_GR(opt.nx, opt.nlab, statename=opt.statename, ncut0=opt.ncut0, ncut1=opt.ncut1, nsig0=opt.nsig0, nsig1=opt.nsig1, fmin=opt.fmin, fmax=opt.fmax, m1=opt.m1, m2=opt.m2, Tfrac=opt.Tfrac, beta=opt.beta, sig=opt.sig)
    
    #plsim.plot_sim(opt.statename+'.npy', './figures/amp_state.png', r'\tilde{A}(f)', fmin=opt.fmin, fmax=opt.fmax, m1=opt.m1, m2=opt.m2, Tfrac=opt.Tfrac, beta=opt.beta, sig=opt.sig)
