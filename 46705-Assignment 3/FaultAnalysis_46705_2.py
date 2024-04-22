"""
46705 - Power Grid Analysis
This file contains the definitions of the functions needed to
carry out Fault Analysis calculations in python.
"""

import numpy as np
import math

# 1. the FaultAnalysis() function
def FaultAnalysis(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,fault_type,Zf,Vf):
    ''' String here with purpose '''
    # calculate sequence fault currents
    Iseq = Calculate_Sequence_Fault_Currents(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,fault_type,Zf,Vf)
    # calculate sequence fault voltages
    Vseq_mat = Calculate_Sequence_Fault_Voltages(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,Vf,Iseq)
    # convert sequence currents to phase (fault) currents
    Iph = Convert_Sequence2Phase_Currents(Iseq)
    # convert sequence voltages to phase line-to-ground (fault) voltages
    Vph_mat = Convert_Sequence2Phase_Voltages(Vseq_mat)    
    return Iph, Vph_mat

# 1.1. the Calculate_Sequence_Fault_Currents() function
def Calculate_Sequence_Fault_Currents(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,fault_type,Zf,Vf):
#  fault_type: 0 = 3-phase balanced fault; 1 = Single Line-to-Ground fault;
#              2 = Line-to-Line fault;     3 = Double Line-to-Ground fault.    
    # Iseq current array: 
    # Iseq[0] = zero-sequence; Iseq[1] = positive-sequence; Iseq[2] = negative-sequence
    Iseq = np.zeros(3,dtype=complex)
    fb = bus_to_ind[fault_bus]
    Znn0=Zbus0[fb][fb]
    Znn1=Zbus1[fb][fb]
    Znn2=Zbus2[fb][fb]
    
    if fault_type == 0:
        Iseq[0]==Iseq[2]==0 #Zero and negative are zero in 3-phase balanced fault
        Iseq[1]=Vf/Znn1 #Znn1 is the bus impedance in positive sequence network
        
    elif fault_type == 1:
        Iseq[0]=Iseq[1]=Iseq[2]=Vf/(Znn0+Znn1+Znn2+3*Zf) #Same sequence currents
        
    elif fault_type == 2:
        Iseq[0]=0
        Iseq[1]=Iseq[2]=Vf/(Znn1+Znn2+Zf) #Same currents in positive and negative sequence
        
    elif fault_type == 3:
        b=(Znn2*(Znn0+3*Zf))/(Znn2+Znn0+3*Zf)
        Iseq[1]=Vf/(Znn1+b)
        Iseq[0]=Iseq[2]=-Iseq[1]
        
    else:
        print('Unknown Fault Type')
    
    return Iseq

# 1.2 the Calculate_Sequence_Fault_Voltages() function
def Calculate_Sequence_Fault_Voltages(Zbus0,Zbus1,Zbus2,bus_to_ind,fault_bus,Vf,Iseq):
    Vseq_mat = np.zeros((len(Zbus0),len(Iseq)),dtype=complex)
    fb = bus_to_ind[fault_bus]
    for k in range (len(Zbus0)):
        Vk0=0-Zbus0[k][fb]*Iseq[0]
        Vk1=Vf-Zbus1[k][fb]*Iseq[1]
        Vk2=0-Zbus2[k][fb]*Iseq[2]
        Vseq_mat[k]=Vk0,Vk1,Vk2
        
    #print(Vseq_mat)
        
    return Vseq_mat

# 1.3. the Convert_Sequence2Phase_Currents() function
def Convert_Sequence2Phase_Currents(Iseq):
    a=-0.5+math.sqrt(3)/2j
    A=[[1,1,1],[1,a**2,a],[1,a,a**2]]
    Iph=np.dot(A,Iseq)
    # enter your code here
    return Iph

# 1.4 the Convert_Sequence2Phase_Voltages() function
def Convert_Sequence2Phase_Voltages(Vseq_mat):
    a=-0.5+math.sqrt(3)/2j
    A=[[1,1,1],[1,a**2,a],[1,a,a**2]]
    N=len(Vseq_mat)
    Vph_mat = np.zeros((N,3),dtype=complex)
    for i in range(N):
        Vph_mat[i]=np.transpose(np.dot(A,np.transpose(Vseq_mat[i])))
        
    return Vph_mat

# ####################################################
# #  Displaying the results in the terminal window   #
# ####################################################
# 2. the DisplayFaultAnalysisResults() function
def DisplayFaultAnalysisResults(Iph,Vph_mat,fault_bus,fault_type,Zf,Vf):
    print('==============================================================')
    print('|                  Fault Analysis Results                    |')
    print('==============================================================')
    # enter your code here
    print('==============================================================')  
    return