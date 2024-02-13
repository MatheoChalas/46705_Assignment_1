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
    bus_data,load_data,gen_data,line_data, tran_data,mva_base, bus_to_ind, ind_to_bus = \
    read_network_data_from_file(file_name)

    MVA_base = mva_base   #OBS...... should be a part of the network data....
    N = len(bus_data)
    Ybus = np.zeros((N,N),dtype=complex)
    
    for line in line_data:
        bus_fr, bus_to, id_, R,X,B_2 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Z_se = R + 1j*X; Y_se = 1/Z_se
        Y_sh_2 = 1j*B_2/2
        
        #Update the matrix:
        Ybus[ind_fr,ind_fr]+= Y_se + Y_sh_2
        Ybus[ind_to,ind_to]+= Y_se + Y_sh_2
        Ybus[ind_fr,ind_to]+= -Y_se
        Ybus[ind_to,ind_fr]+= -Y_se
        
        
    #Get transformer data as well...
    bus_kv = []
    buscode = []
    bus_labels = []
    for line in bus_data:
        b_nr, label, kv, code = line
        buscode.append(code)
        bus_labels.append(label)
        bus_kv.append(kv)
        
    buscode = np.array(buscode)
    bus_kv = np.array(bus_kv)


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

#bus-branch matrices
    N_branches = len(line_data) + len(tran_data)
    br_f = -np.ones(N_branches,dtype=int)
    br_t = -np.ones(N_branches,dtype=int)
    
    Y_fr = np.zeros((N_branches,N),dtype=complex)
    Y_to = np.zeros((N_branches,N),dtype=complex)
        
    for line,i in zip(line_data,range(len(line_data))):
        bus_fr, bus_to, id_, R,X,B_2 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Z_se = R + 1j*X; Y_se = 1/Z_se
        Y_sh_2 = 1j*B_2/2
        # update the entries
        Y_fr[i,ind_fr] =  Y_se + Y_sh_2       
        Y_fr[i,ind_to] = -Y_se
        Y_to[i,ind_to] =  Y_se + Y_sh_2       
        Y_to[i,ind_fr] = -Y_se
        br_f[i] = ind_fr
        br_t[i] = ind_to
    
    for line,i in zip(tran_data,range(len(line_data),N_branches)):
        bus_fr, bus_to, id_, R,X,n,ang1 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Zeq = R+1j*X; Yeq = 1/Zeq
        c = n*np.exp(1j*ang1/180*np.pi)
        ### adminttance matrix
        Yps_mat = np.zeros((2,2),dtype=complex)
        Yps_mat[0,0] = Yeq/np.abs(c)**2
        Yps_mat[0,1] = -Yeq/c.conj()
        Yps_mat[1,0] = -Yeq/c
        Yps_mat[1,1] = Yeq
        # indices
        ind_ = np.array([ind_fr,ind_to])
        #update
        Ybus[np.ix_(ind_,ind_)] += Yps_mat
        br_f[i] = ind_fr
        br_t[i] = ind_to
        # update the entries
        Y_fr[i,ind_fr] =  Yps_mat[0,0]      
        Y_fr[i,ind_to] =  Yps_mat[0,1]
        Y_to[i,ind_to] =  Yps_mat[1,1]       
        Y_to[i,ind_fr] =  Yps_mat[1,0]
        
  return Ybus, Sbus, V0,pv_index,pq_index

print(LoadNetworkData(filename))

