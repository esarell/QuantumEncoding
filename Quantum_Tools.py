from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, Aer, IBMQ,transpile
#from scipy.interpolate import approximate_taylor_polynomial
#from qiskit.circuit.library import RGQFTMultiplier, DraperQFTAdder, ExactReciprocal
from qiskit.circuit.library.basis_change import QFT as QFT_pre
#from qiskit.extensions import HamiltonianGate
#from qiskit.circuit.library.standard_gates import PhaseGate, RYGate, CSwapGate
#from qiskit.circuit.library.arithmetic.integer_comparator import IntegerComparator
from qiskit.circuit.library.standard_gates import XGate, PhaseGate
from qiskit.quantum_info import Pauli
from qiskit.circuit.library.boolean_logic import OR
#from qiskit.quantum_info import random_unitary
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Unroller
import itertools as it
import numpy as np
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
    print("QFT")
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
    print("n1:",n1)
    print("binary:",binary)
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
                print("check")

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
    Args:
    digit: input number (float)
    n: total number of bits (int)
    nint: number of integer bits. Default to lowest required (int)
    phase: handle negative numbers
    Return:
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
    #print(compliment)
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
                print("append")
        circ.x(qreg[0])
        print("or here")

    else:
        '''Currently this bit of code runs every time'''
        #Adds one to the current value in the qlab
        bin_one = ''.join(np.zeros(n-1).astype(int).astype(str))+'1'
        inc_gate = QFTBinaryAdd(circ, qreg, bin_one, wrap=True, QFT_on=QFT_on, iQFT_on=iQFT_on)
        circ.append(inc_gate, qreg)
        #print("here")

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'

    return circ

def integer_compare(circ, qreg, qtarg, qans ,value, geq=True, wrap=True, inverse=False, uncomp=True, label='P'):
    
    n = len(qreg)
    size_anc =len(qans)

    if wrap:
        qreg = QuantumRegister(n, 'q_reg')
        qans = QuantumRegister(size_anc, 'q_ans')
        qtarg = QuantumRegister(1,'qtarg')
        circ = QuantumCircuit(qreg,qtarg, qans)
    if value<=0.:
        #geq is currently set to False
        circ.x(qtarg)
        if False:
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
                        #AND GATE
                        circ.ccx(qreg[i], qans[i - 1], qans[i]);
                #Once we have reached the last bit we do the same checks but now but the carry bit
                #In the target register instead of the ancilary register
                else:
                    if twos[i]=='1':
                        circ.compose(OR(2), [qreg[i], qans[i-1], qtarg[0]], inplace=True);
                    else:
                        circ.ccx(qreg[i], qans[i - 1], qtarg[0]);
            if not geq:
               circ.x(qtarg);
            
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
        
            if not geq:
                circ.x(qtarg);
    else:
        if not geq:
            circ.x(qtarg);
    
    if wrap:
        circ = circ.to_gate()
        circ.label = label


    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'
    
    return circ

def QFTMultiply(circ, qreg1, qreg2, qreg3, A=1., wrap=False, inverse=False, nint1=None, nint2=None, nint3=None, phase=False, label='Mult', QFT_on=True, iQFT_on=True):
    """
    |qreg1>|qreg2>|qreg3> -> |qreg1>|qreg2>|qreg1 x qreg2>
    """
    n1 = len(qreg1)
    n2 = len(qreg2)
    n3 = len(qreg3)
    print("n1",n1)
    print("n2",n2)
    print("n3",n3)
    if n3!=n1+n2 and nint3 == None:
        raise ValueError('Output register should be the combined length of both input registers if no integer bit length is specified.')

    if nint1==None:
        nint1=n1
    if nint2==None:
        nint2=n2
    if nint3==None:
        nint3=n3
    
    nshift = (nint1 + nint2)-nint3

    if phase:
        nshift+=1

    if inverse:
        wrap = True

    if wrap:
        qreg1 = QuantumRegister(n1, 'q_reg1')
        qreg2 = QuantumRegister(n2, 'q_reg2')
        qreg3 = QuantumRegister(n3, 'q_reg3')
        circ = QuantumCircuit(qreg1, qreg2, qreg3)

    if QFT_on:
        circ.append(QFT(circ, qreg3, do_swaps=False, wrap=True), qreg3[:])

    for j in np.arange(1, n1 + 1):
        for i in np.arange(1, n2 + 1):
            for k in np.arange(1, n3 + 1):
                lam = A*(2 * np.pi) / (2. ** (i + j + k - n3 - nshift))
                if lam%(2*np.pi)==0.:
                    continue
                circ.append(PhaseGate(lam).control(2),[qreg1[n1 - j], qreg2[n2 - i], qreg3[k - 1]])

    if iQFT_on:
        circ.append(QFT(circ, qreg3, do_swaps=False, wrap=True, inverse=True), qreg3[:])

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'\dag'

    return circ
def QFTAddition_(circ, qreg1, qreg2, wrap=False, inverse=False, label='Add', QFT_on=True, iQFT_on=True, pm=1, phase=False):
    r"""
    |qreg1>|qreg2> -> |qreg1>|qreg2 + qreg1>
    """
    n1 = len(qreg1)
    n2 = len(qreg2)

    if inverse:
        wrap = True

    if wrap:
        qreg1 = QuantumRegister(n1, 'q_reg1')
        qreg2 = QuantumRegister(n2, 'q_reg2')
        circ = QuantumCircuit(qreg1, qreg2)

    if QFT_on:
        circ.append(QFT(circ, qreg2, do_swaps=False, wrap=True), qreg2[:])

    jend = n1
    if phase and n1!=n2:
        circ.append(twos_compliment(circ, qreg1[:-1], wrap=True).control(1), [qreg1[-1], *qreg1[:-1]]);
        jend -= 1

    for j in np.arange(0,jend):
        for l in np.arange(0,n2):
            lam = pm*2*np.pi*2.**(j-l-1)
            if lam%(2*np.pi)==0:
                continue
            circ.cp(lam, qreg1[j], qreg2[l])
            if phase and n1!=n2 and lam%np.pi!=0:
                circ.append(PhaseGate(-2.*lam).control(2), [qreg1[-1], qreg1[j], qreg2[l]]);

    if phase and n1!=n2:
        circ.append(twos_compliment(circ, qreg1[:-1], wrap=True).control(1), [qreg1[-1], *qreg1[:-1]]);

    if iQFT_on:
        circ.append(QFT(circ, qreg2, do_swaps=False, wrap=True, inverse=True), qreg2[:])

    if wrap:
        circ = circ.to_gate()
        circ.label = label

    if inverse:
        circ = circ.inverse()
        circ.label = label+'†'

    return circ


if __name__ == "__main__":
    n=4
    qreg = QuantumRegister(n, 'q_reg')
    qans = QuantumRegister(n, 'q_ans')
    qtarg = QuantumRegister(1,'qtarg')
    circ = QuantumCircuit(qreg,qtarg, qans)
    integer_circ = integer_compare(circ,qreg,qtarg,qans,1,inverse=False)
    circ=circ.compose(integer_circ,[*qreg,*qtarg,*qans])
    circ.decompose().draw("mpl")
    plt.show()