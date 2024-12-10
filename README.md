# QuantumEncoding
This repository contains two projects, Reduced Amplitude Encoding using fixed rotations (Fixed_Rotations_GSP) and 
Grover-Algorithm using Linear Piecewise (GroverRudolph_LPW).

## Set Up
```
pip -r requirements.txt
```

# Reduced Amplitude Encoding using Fixed Rotations
Attempting to reduce the amount of Rotational Y gates needed to encode a distribution using the Grover-Rudolph algorithm.
Marin-Sanchez, et al., 2023, showed that you can reduce a circuit to a single rotation for later on levels:
https://link.aps.org/doi/10.1103/PhysRevResearch.5.033114
<br />
We adapt this to have reduced rotations but with additional conditons, so that we have fixed rotations on specific
sections of our distribution.

## Analytical Solution (Inspiral)
To test for an analytical solution we looked at encoding the Inspiral of a Gravitational Wave. We can approximate the amplitude
of this function to be  $`A(f ) ∝ f^{−\frac{7}{6}}`$ 
### inspiral_cal_k.py
This python file calculates the $`\eta `$ values for each level of the algoirthm, and then calculates a correponding $`k `$ value.
This tells us the level we could switch to a fixed rotation for that section.
<br />
To change any of the paramaters for this code simply change the main() function.
It is currently set to 9 qubits with an error rate of 0.9999. 
### fixed_rotation_circ.py
This program creates a quantum circuit which encodes the inspiral amplitudes.
It then calculates the fidelity of the quantum state to the actual amplitude values.
## Numerical Solution (Waveform)
To test this further we encoded a gravitational wave which doesn't have an analytical solution. 
### Numerical Data
This is generated using the LinearSplineAprrox.py file.
<br />
#### waveform_data.pkl
This is a pickle file which contains all the information on the waveform for 512 datapoints. 
It has a frequency range of 40-350Hz.
#### waveform_amps.pkl
This is a pickle file containing the amplitudes of the waveform for 512 datapoints. 

## waveform_fixed_rotation_cal.py
This code performs the second devirative of the log of the function squared.
For the whole waveform this is done numerically done using points which are generated in LinearSpineApprox.py
Once this $`\eta `$ value is calculated we determine the value of $`k `$ which tells us the level at which we switch to a fixed rotational gate

## wave_fixed_circ.py
This program creates a quantum circuit which encodes the waveform amplitudes.
It then calculates the fidelity of the quantum state to the actual amplitude values.

# Grover-Algorithm using Linear Piecewise
Work in progress.
## Using qiskit_tools.py
This program is based off the paper by Hayes, et al.
https://arxiv.org/abs/2306.11073
And makes uses of the qiskit_tools.py file which can be found:
https://github.com/Fergus-Hayes/qiskit_tools
## How to Run
Start by installing the requirements.txt 
Then you should be able to run the GroverRudolph.py file.
It has a generic set up currently and needs to be genralised, it currently is working for 4 qubits. 
It uses 9 qubits for coeffeicents.
3 for labels
9 for ancillary qubits
1 for the target bit
Which is a total of 26 qubits however the target bit could be included in the ancillary bits to reduce the amoutn of qubits needed
Currently this is using the data that Ashwin generated which can be seen in the GenBounds.py file
## Key files
### Quantum_Tools.py
contains the key functions from Fergus' qiskit_tools.py file with comments and it is what the other files call.
### Linear_Piecewise.py 
contains all the functions for the linear piecewise function, including the label gate, the loading and unloading of coefficents. 
### GenBounds.py 
is a modifided version of the jyupter notebook provided by Ashwin, it generates the data points using the function and the second derivative. It uses the amount of qubits we have to decide how many bounds we should have

### GroverRudolph.py
