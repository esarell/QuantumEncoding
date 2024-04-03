# QuantumEncoding
To put all my Grover-Rudolph algorithm code in. Will probably end up being the main repository for all amplitude encoding things
Test
## Using qiskit_tools.py
This program is based off the paper by Hayes, et al.
https://arxiv.org/abs/2306.11073
And makes uses of the qiskit_tools.py file which can be found:
https://github.com/Fergus-Hayes/qiskit_tools
# How to Run
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