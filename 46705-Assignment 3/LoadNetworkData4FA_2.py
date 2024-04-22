# -*- coding: utf-8 -*-
import numpy as np
import ReadNetworkData4FA as rd4fa

def LoadNetworkData4FA(filename):
    global Ybus,Sbus,V0,buscode,pq_index,pv_index,Y_fr,Y_to,br_f,br_t,br_Y,S_LD, \
           ind_to_bus,bus_to_ind,MVA_base,bus_labels,Ybus0,Ybus2,Zbus0,Zbus1,Zbus2          
    # read in the data from the file...
    bus_data,load_data,gen_data,line_data,tran_data,mva_base,bus_to_ind,ind_to_bus = \
    rd4fa.read_network_data_4fa_from_file(filename)

    ############################################################################################## 
    # Construct the Ybus (positive-sequence), Ybus0 (zero-sequence), and Ybus2 (negative-sequence)
    # matrices from elements in the line_data and trans_data
    # Keep/modify code from the Python power flow program as needed
    ##########################################################################  
    N = len(bus_data) # Number of buses
    Ybus = np.zeros((N,N),dtype=complex) #Positive sequence
    N_branches = len(line_data) + len(tran_data)
    #bus-branch matrices
    N_branches = len(line_data) + len(tran_data)
    br_f = -np.ones(N_branches,dtype=int)
    br_t = -np.ones(N_branches,dtype=int)
    
#Ybus - Positive sequence
    #Get data from lines
    for line in line_data:
        bus_fr, bus_to, id_, R,X,B_2,X2, X0 = line #unpack the lines
        #Computation of the index and the impedance and admittance of each lines
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Z_se = 1j*X; Y_se = 1/Z_se #Resistance and shunt are neglected 
        
        #Update the matrix:
        Ybus[ind_fr,ind_fr]+= Y_se 
        Ybus[ind_to,ind_to]+= Y_se 
        Ybus[ind_fr,ind_to]+= -Y_se
        Ybus[ind_to,ind_fr]+= -Y_se
        
    #Get data from transformers   
    for line,i in zip(tran_data,range(len(line_data),N_branches)):
        bus_fr, bus_to, id_, R,X,n,ang1,fr_co, to_co, X2, X0 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Zeq = 1j*X; Yeq = 1/Zeq #We neglect the resistance
        #c = n*np.exp(1j*ang1/180*np.pi) we neglect the sifhtangle
        Yps_mat = np.zeros((2,2),dtype=complex)
        ### adminttance matrix
        Yps_mat[0,0] = Yeq
        Yps_mat[0,1] = -Yeq
        Yps_mat[1,0] = -Yeq
        Yps_mat[1,1] = Yeq
        # indices
        ind_ = np.array([ind_fr,ind_to])
        #update
        Ybus[np.ix_(ind_,ind_)] += Yps_mat
        br_f[i] = ind_fr
        br_t[i] = ind_to
        # update the entries
        #Y_fr[i,ind_fr] =  Yps_mat[0,0]      
        #Y_fr[i,ind_to] =  Yps_mat[0,1]
        #Y_to[i,ind_to] =  Yps_mat[1,1]       
        #Y_to[i,ind_fr] =  Yps_mat[1,0]
        

        #br_Ymat[i]=Yps_mat

    #Get data from generators
    for line in gen_data:
        bus_nr, MVA_size, p_gen,X, X2, X0, Xn, GRND = line #new     
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Zeq = 1j*X; Yeq = 1/Zeq #We neglect the resistance
        Yps_mat = np.zeros((2,2),dtype=complex)
        ### adminttance matrix
        Yps_mat[0,0] = Yeq
        Yps_mat[0,1] = -Yeq
        Yps_mat[1,0] = -Yeq
        Yps_mat[1,1] = Yeq
        # indices
        ind_ = np.array([ind_fr,ind_to])
        #update
        Ybus[np.ix_(ind_,ind_)] += Yps_mat
        br_f[i] = ind_fr
        br_t[i] = ind_to

#Ybus2 - Negative sequence      
    Ybus2 = np.zeros((N,N),dtype=complex) #Negative sequence   
    #Get data from lines
    for line in line_data:
        bus_fr, bus_to, id_, R,X,B_2,X2, X0 = line #unpack the lines
        #Computation of the index and the impedance and admittance of each lines
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Z_se = 1j*X2; Y_se = 1/Z_se #Resistance and shunt are neglected. We consider X2
        
        #Update the matrix:
        Ybus2[ind_fr,ind_fr]+= Y_se 
        Ybus2[ind_to,ind_to]+= Y_se 
        Ybus2[ind_fr,ind_to]+= -Y_se
        Ybus2[ind_to,ind_fr]+= -Y_se
        
    #Get data from transformers   
    for line,i in zip(tran_data,range(len(line_data),N_branches)):
        bus_fr, bus_to, id_, R,X,n,ang1,fr_co, to_co, X2, X0 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Zeq = 1j*X2; Yeq = 1/Zeq #We neglect the resistance. We consider X2
        #c = n*np.exp(1j*ang1/180*np.pi) we neglect the sifhtangle
        Yps_mat = np.zeros((2,2),dtype=complex)
        ### adminttance matrix
        Yps_mat[0,0] = Yeq
        Yps_mat[0,1] = -Yeq
        Yps_mat[1,0] = -Yeq
        Yps_mat[1,1] = Yeq
        # indices
        ind_ = np.array([ind_fr,ind_to])
        #update
        Ybus2[np.ix_(ind_,ind_)] += Yps_mat
        br_f[i] = ind_fr
        br_t[i] = ind_to
        # update the entries
        #Y_fr[i,ind_fr] =  Yps_mat[0,0]      
        #Y_fr[i,ind_to] =  Yps_mat[0,1]
        #Y_to[i,ind_to] =  Yps_mat[1,1]       
        #Y_to[i,ind_fr] =  Yps_mat[1,0]
        

        #br_Ymat[i]=Yps_mat

    #Get data from generators
    for line in gen_data:
        bus_nr, MVA_size, p_gen,X, X2, X0, Xn, GRND = line #new     
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Zeq = 1j*X2; Yeq = 1/Zeq #We neglect the resistance. We consider X2
        Yps_mat = np.zeros((2,2),dtype=complex)
        ### adminttance matrix
        Yps_mat[0,0] = Yeq
        Yps_mat[0,1] = -Yeq
        Yps_mat[1,0] = -Yeq
        Yps_mat[1,1] = Yeq
        # indices
        ind_ = np.array([ind_fr,ind_to])
        #update
        Ybus2[np.ix_(ind_,ind_)] += Yps_mat
        br_f[i] = ind_fr
        br_t[i] = ind_to
    
 #Ybus0 - Zero sequence    
    Ybus0 = np.zeros((N,N),dtype=complex) #Zero sequence
    #Get data from lines
    for line in line_data:
        bus_fr, bus_to, id_, R,X,B_2,X2, X0 = line #unpack the lines
        #Computation of the index and the impedance and admittance of each lines
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to] 
        Z_se = 1j*X0; Y_se = 1/Z_se #Resistance and shunt are neglected. We consider X0
        
        #Update the matrix:
        Ybus0[ind_fr,ind_fr]+= Y_se 
        Ybus0[ind_to,ind_to]+= Y_se 
        Ybus0[ind_fr,ind_to]+= -Y_se
        Ybus0[ind_to,ind_fr]+= -Y_se
        
    #Get data from generators
    for line in gen_data:
        bus_nr, MVA_size, p_gen,X, X2, X0, Xn, GRND = line #new     
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to]
        if GRND==1:
            Zeq = 1j*(X0+3*Xn); Yeq = 1/Zeq #We neglect the resistance. We consider X0 and Xn
        else:
            Zeq = 1j*(X0); Yeq = 1/Zeq #We neglect the resistance. We consider X0 (Xn generates an open circuit)
        Yps_mat = np.zeros((2,2),dtype=complex)
        ### adminttance matrix
        Yps_mat[0,0] = Yeq
        Yps_mat[0,1] = -Yeq
        Yps_mat[1,0] = -Yeq
        Yps_mat[1,1] = Yeq
        # indices
        ind_ = np.array([ind_fr,ind_to])
        #update
        Ybus0[np.ix_(ind_,ind_)] += Yps_mat
        br_f[i] = ind_fr
        br_t[i] = ind_to
    
    #Get data from transformers   
    for line,i in zip(tran_data,range(len(line_data),N_branches)):
        bus_fr, bus_to, id_, R,X,n,ang1,fr_co, to_co, X2, X0 = line #unpack
        ind_fr = bus_to_ind[bus_fr]    
        ind_to = bus_to_ind[bus_to]
        if fr_co==2 and to_co==1:
            Zeq = 1j*X0; Yeq = 1/Zeq #We neglect the resistance. We consider X0
        elif fr_co==2 and to_co==2:
            Zeq = 1j*X0; Yeq = 1/Zeq #We neglect the resistance. We consider X0
        elif fr_co==2 and to_co==3:
            Zeq = 1j*X0; Yeq = 1/Zeq #We neglect the resistance. We consider X0
        elif fr_co==1 and to_co==3:
            Zeq = 1j*X0; Yeq = 1/Zeq #We neglect the resistance. We consider X0
        elif fr_co==3 and to_co==3:
            Zeq =0; Yeq = 0 #We neglect the resistance. We consider X0
            

        #c = n*np.exp(1j*ang1/180*np.pi) we neglect the sifhtangle
        Yps_mat = np.zeros((2,2),dtype=complex)
        ### adminttance matrix
        Yps_mat[0,0] = Yeq
        Yps_mat[0,1] = -Yeq
        Yps_mat[1,0] = -Yeq
        Yps_mat[1,1] = Yeq
        # indices
        ind_ = np.array([ind_fr,ind_to])
        #update
        Ybus0[np.ix_(ind_,ind_)] += Yps_mat
        br_f[i] = ind_fr
        br_t[i] = ind_to
        # update the entries
        #Y_fr[i,ind_fr] =  Yps_mat[0,0]      
        #Y_fr[i,ind_to] =  Yps_mat[0,1]
        #Y_to[i,ind_to] =  Yps_mat[1,1]       
        #Y_to[i,ind_fr] =  Yps_mat[1,0]
        

        #br_Ymat[i]=Yps_mat
    Zbus0 = np.linalg.inv(Ybus0)
    Zbus1 = np.linalg.inv(Ybus)
    Zbus2 = np.linalg.inv(Ybus2)
    
    return   