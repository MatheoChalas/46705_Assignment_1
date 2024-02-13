import numpy as np
import ReadNetworkData as rd


"""Format of variables"""

"""
bus_data: rows contain [bus_nr, bus_label, Kv_level, code]
load_data: rows contain [bus_nr, P_LD_MW, Q_LD_MW] (loads in MW)
gen_data: rows contain [bus_nr, MVA_size, P_Gen] (in MVA and MW)
line_data: rows contain [fr_bus, to_bus, ID_label, R, X, B/2] (pu system base)
tran_data: contains [fr_bus, to_bus, ID_label, R_eq, X_eq, n_pu, ang_deg]
(R_eq and X_eq in pu on system base, n_pu: the pu turns ratio, ang_deg: is the phase shift)
mva_base: the system MVA base
bus_to_ind: mapping from bus numbers to the corresponding indices in the bus matrices and arrays
ind_to_bus: containing the mapping from the indices to the busses (the opposite of the above)
"""
filename = 'TestSystem.txt'

def LoadNetworkData(filename):
  global Ybus, Y_from, Y_to, br_f, br_t, buscode, bus_labels, S_LDMVA_base
  
  #read in the data from the file...
  bus_data, load_data, gen_data, line_data, tran_data, mva_base, bus_to_ind, ind_to_bus = rd.read_network_data_from_file(filename)
  
  N=len(line_data[0])
  Ybus = np.zeros((N,N))
  Y=np.array([1/(line_data[k,3]+1.J*line_data[k,4]) for k in range(N)])
  B=np.array([(1.J*line_data[k,5])/2 for k in range(N)])
  
  for k in range(N):
    for i in range(N):
      if k==i:
        Ybus[k,i]=sum(Y,B,-Y[k]-B[k])
      if k<i:
        Ybus[k,i]=-Y[i-1]
      else :
        Ybus[k,i]=Ybus[i,k]

 Sbus= np.zeros(N, dtype=complex)
 SLD = np.zeros(N, dtype=complex)
 
 for line in load_data:
     bus_nr, PLD, QLD = line
     ind_nr = bus_to_ind[bus_nr]
     SLD =(PLD+1j*QLD)/MVA_base
     Sbus[ind_nr] += -SLD 
     
 for line in gen_data:
     bus_nr, MVA_size, p_gen = line
     ind_nr = bus_to_ind[bus_nr]
     SLD =(p_gen)/MVA_base
     Sbus[ind_nr] += SLD 

V0 = np.ones(N,dtype=np.complex)
buscode = np.array([3, 2, 1])
pq_index = np.where(buscode == 1)[0]
pv_index = np.where(buscode == 2)[0]
ref = np.where(buscode == 3)[0]
  return Ybus, Sbus, SLD

print(LoadNetworkData(filename))

