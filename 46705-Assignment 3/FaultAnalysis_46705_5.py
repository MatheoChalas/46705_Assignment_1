"""
46705 - Power Grid Analysis
This file contains the definitions of the functions needed to
carry out Fault Analysis calculations in python.
"""

import numpy as np
import math
import cmath

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
        Iseq[0]==Iseq[2]==0 
        Iseq[1]=Vf/Znn1 
        
    elif fault_type == 1:
        Iseq[0]=Iseq[1]=Iseq[2]=Vf/(Znn0+Znn1+Znn2+3*Zf) 
        
    elif fault_type == 2:
        Iseq[0]=0
        Iseq[1]=Vf/(Znn1+Znn2+Zf)
        Iseq[2]=-Iseq[1]
        
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
        Vk0=-Zbus0[k][fb]*Iseq[0]
        Vk1=Vf-Zbus1[k][fb]*Iseq[1]
        Vk2=-Zbus2[k][fb]*Iseq[2]
        Vseq_mat[k]=Vk0,Vk1,Vk2
       
    return Vseq_mat

# 1.3. the Convert_Sequence2Phase_Currents() function
def Convert_Sequence2Phase_Currents(Iseq):
    a=-0.5+ 1j*(math.sqrt(3)/2)
    A=[[1,1,1],[1,a**2,a],[1,a,a**2]]
    Iph=np.dot(A,Iseq)
    # enter your code here
    return Iph

# 1.4 the Convert_Sequence2Phase_Voltages() function
def Convert_Sequence2Phase_Voltages(Vseq_mat):
    a=-0.5+ 1j*(math.sqrt(3)/2)
    A=[[1,1,1],[1,a**2,a],[1,a,a**2]]
    N=len(Vseq_mat)
    Vph_mat = np.zeros((N,3),dtype=complex)
    for i in range(N):
        Vph_mat[i]=np.dot(A,(Vseq_mat[i]))
    return Vph_mat

# ####################################################
# #  Displaying the results in the terminal window   #
# ####################################################
# 2. the DisplayFaultAnalysisResults() function
def DisplayFaultAnalysisResults(Iph,Vph_mat,fault_bus,fault_type,Zf,Vf):
    fault_name=""
    if fault_type == 0:
        fault_name="| 3-phase balanced fault at Bus %d." %(fault_bus)
    elif fault_type == 1:
        fault_name="| Single Line-to-Ground fault at Bus %d, phase a." %(fault_bus)
    elif fault_type == 2:
        fault_name="| Line-to-Line fault at Bus %d, between phase b and phase c." %(fault_bus)
    elif fault_type == 3:
        fault_name="| Double Line-to-Ground fault at Bus %d, phase b and c." %(fault_bus)
   
    Ia=Iph[0]
    Ib=Iph[1]
    Ic=Iph[2]    
    
    data=[]
    for k in range(len(Vph_mat)):
        line=[k+1,abs(Vph_mat[k,0]),cmath.phase(Vph_mat[k,0])*180/math.pi,abs(Vph_mat[k,1]),cmath.phase(Vph_mat[k,1])*180/math.pi,abs(Vph_mat[k,2]),cmath.phase(Vph_mat[k,2])*180/math.pi]
        data.append(line)
    
    print('==============================================================')
    print('|                  Fault Analysis Results                    |')
    print('==============================================================')
    # enter your code here
    print(fault_name)
    print('| Prefault Voltage: Vf = %.3f  (pu)                         |'%(Vf))
    print('| Fault Impedance:  Zf = %.3f  (pu)                         |'%(Zf))  
    print('==============================================================')
    print('| Phase Currents --------------------------------------------|')
    print('| ----------------                                           |')
    print('|     ---- Phase a ----| ---- Phase b ----| ---- Phase c ----|')
    print('|     -----------------|------------------|------------------|')
    print('|      Mag(pu) Ang(deg)|  Mag(pu) Ang(deg)|  Mag(pu) Ang(deg)|')
    print('|       %.3f   %.2f |   %.3f  %.2f  |   %.3f  %.2f   |'%(abs(Ia),cmath.phase(Ia)*180/math.pi,abs(Ib),cmath.phase(Ib)*180/math.pi,abs(Ic),cmath.phase(Ic)*180/math.pi))
    
    print("==============================================================")
    print("| Phase Line-to-Ground Voltages -----------------------------|")
    print("==============================================================")
    print("| -----------------------------                              |") 
    print("|   | ---- Phase a ----| ---- Phase b ----| ---- Phase c ----|") 
    print("|Bus|------------------|------------------|------------------|") 
    print("|   |  Mag(pu) Ang(deg)|  Mag(pu) Ang(deg)|  Mag(pu) Ang(deg)|") 
    print("|---| -------- ------- | -------- ------- | -------- ------- |")
    print("| %d |   %.3f     %.2f |  %.3f  %.2f  |   %.3f  %.2f  |"%(data[0][0],data[0][1],data[0][2],data[0][3],data[0][4],data[0][5],data[0][6])) 
    print("| %d |   %.3f     %.2f |  %.3f  %.2f  |   %.3f  %.2f  |"%(data[1][0],data[1][1],data[1][2],data[1][3],data[1][4],data[1][5],data[1][6]))
    print("| %d |   %.3f     %.2f |  %.3f   %.2f    |   %.3f  %.2f    |"%(data[2][0],data[2][1],data[2][2],data[2][3],data[2][4],data[2][5],data[2][6]))
    print("| %d |   %.3f     %.2f |  %.3f  %.2f  |   %.3f  %.2f  |"%(data[3][0],data[3][1],data[3][2],data[3][3],data[3][4],data[3][5],data[3][6]))
    print("| %d |   %.3f     %.2f |  %.3f  %.2f  |   %.3f  %.2f  |"%(data[4][0],data[4][1],data[4][2],data[4][3],data[4][4],data[4][5],data[4][6]))
    print("==============================================================") 


    
  
    return