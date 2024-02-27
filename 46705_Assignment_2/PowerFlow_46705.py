"""
46705 - Power Grid Analysis
This file contains the definitions of the functions needed to
carry out Power Flow calculations in python.

How to carry out Power Flow in a new *.py file? 
See the example in table 1 in the assignment text
"""

import numpy as np
from tabulate import tabulate

# 1. the PowerFlowNewton() function
def PowerFlowNewton(Ybus,Sbus,V0,pv_index,pq_index,max_iter,err_tol):
    ''' String here with purpose '''
    # implement your code here
    success = 0 #Initializationof statusflagand iterationcounter
    n = 0
    V = V0
    print(' iteration maximumP &Q mismatch(pu)')
    print('------------------------------------')
    #Determinemismatchbetweeninitialguessand andspecifiedvalueforP andQ
    F = calculate_F(Ybus,Sbus,V,pv_index,pq_index)
    #Checkif thedesiredtoleranceisreached
    success =CheckTolerance(F,n,err_tol)
    #Startthe Newtoniterationloop
    while (not success) and (n < max_iter):
        n += 1 # Updatecounter
        # ComputederivativesandgeneratetheJacobianmatrix
        J_dS_dVm,J_dS_dTheta = generate_Derivatives(Ybus,V)
        J = generate_Jacobian(J_dS_dVm,J_dS_dTheta,pv_index,pq_index)
        # Computetheupdatestep
        dx = np.linalg.solve(J,F)
        # Updatevoltagesand checkiftoleranceisnow reached
        V = Update_Voltages(dx,V,pv_index,pq_index)
        F = calculate_F(Ybus,Sbus,V,pv_index,pq_index)
        success = CheckTolerance(F,n,err_tol)
    
    if success: #printoutmessageconcerningwetherthe powerflowconvergedornot
        print('The NewtonRapsonPowerFlow Converged in %d iterations!'% (n,))
    else:
        print('No Convergence!!!\n Stopped after %d iterations without solution...'% (n,))
    return V,success,n


# 2. the calculate_F() function
def calculate_F(Ybus,Sbus,V,pv_index,pq_index):
    Delta_P=np.zeros(pv_index+pq_index)
    Delta_Q=np.zeros(pq_index+1)
    
    Delta_S= Sbus- V * (Ybus.dot(V)).conj()
    Delta_P[pv_index]=Delta_S[pv_index].real
    Delta_P[pq_index]=Delta_S[pq_index].real
    Delta_Q[pq_index]=Delta_S[pq_index].imag
    print(Delta_Q)
    F= np.concatenate((Delta_P[pv_index],Delta_P[pq_index],Delta_Q[pq_index]),axis=0)

    return F


# 3. the CheckTolerance() function
def CheckTolerance(F,n,err_tol):
    normF = np.linalg.norm(F,np.inf)
    print("normF= %f"%(normF,),n)
    if normF<err_tol:
        success=1
        
    else :
        success=0
    return success

# 4. the generate_Derivatives() function
def generate_Derivatives(Ybus,V):
    J_ds_dVm=np.diag(V/np.absolute(V)).dot(np.diag((Ybus.dot(V)).conj()))+ np.diag(V).dot(Ybus.dot(np.diag(V/np.absolute(V))).conj())
    J_dS_dTheta = 1j*np.diag(V).dot((np.diag(Ybus.dot(V))-Ybus.dot(np.diag(V))).conj())

    return J_ds_dVm,J_dS_dTheta


# 5. the generate_Jacobian() function
def generate_Jacobian(J_dS_dVm,J_dS_dTheta,pv_index,pq_index):
    #AppendPV andPQ indices for convenience
    pvpq_ind= np.append(pv_index,pq_index)
    
    #Create the sub-matrices
    J_11= np.real(J_dS_dTheta[np.ix_(pvpq_ind,pvpq_ind)])
    J_12= np.real(J_dS_dVm[np.ix_(pvpq_ind,pq_index)])
    J_21= np.imag(J_dS_dTheta[np.ix_(pq_index,pvpq_ind)])
    J_22= np.imag(J_dS_dVm[np.ix_(pq_index,pq_index)])

    #Compute the Jacobian
    J= np.block([[J_11,J_12],[J_21,J_22]])

    return J


# 6. the Update_Voltages() function
def Update_Voltages(dx,V,pv_index,pq_index):
    # Note:differencebetween Python and Matlab when using indices
    N1 = 0; N2 = len(pv_index) # dx[N1:N2]-ang.onthe pvbuses
    N3 = N2; N4 = N3 + len(pq_index)# dx[N3:N4]-ang.onthe pqbuses
    N5 = N4; N6 = N5 + len(pq_index)# dx[N5:N6]-mag.onthe pqbuses

    Theta= np.angle(V);Vm =np.absolute(V)
    if len(pv_index)>0:
        Theta[pv_index]+= dx[N1:N2]
    if len(pq_index)>0:
        Theta[pq_index]+= dx[N3:N4]
        Vm[pq_index]+=dx[N5:N6]
    V = Vm * np.exp(1j*Theta)

    return V



####################################################
#  Displaying the results in the terminal window   #
####################################################
from tabulate import tabulate
def DisplayResults(V,lnd):

    Ybus = lnd.Ybus
    Y_from = lnd.Y_fr
    Y_to = lnd.Y_to
    br_f = lnd.br_f
    br_t = lnd.br_t
    buscode = lnd.buscode
    SLD = lnd.SLD
    ind_to_bus = lnd.ind_to_bus
    bus_to_ind = lnd.bus_to_ind
    MVA_base = lnd.MVA_base
    bus_labels = lnd.bus_labels
    Sbus=lnd.Sbus
    
    Gen_MVA = lnd.Gen_MVA
    
    # Bus results
    S_inj = V*(Ybus.dot(V)).conj()
    bus_results = []
    pq_index = np.where(buscode == 1)[0] # Find indices for all PQ-busses
    pv_index = np.where(buscode == 2)[0] # Find indices for all PV-busses
    ref = np.where(buscode == 3)[0] # Find index for ref bus
    
    #Extraction of the needed values for all the busses 
    k=0
    
    for i in range(len(bus_labels)):
        bus_index = i
        bus_label = bus_labels[i]
        bus_voltage_mag = round(abs(V[bus_index]),3)
        bus_voltage_ang = round(np.angle(V[bus_index], deg=True),2)
        #Initialization of the powers to "-"
        load_P = "-"
        load_Q = "-"
        generation_P = "-"
        generation_Q = "-"
        loading="_"
    
    #Differenciation of the loads and generators
        if buscode[bus_index]==1:
            load_P = -Sbus[bus_index].real
            load_Q = -Sbus[bus_index].imag
            
            
        else:
            generation_P = round(S_inj[bus_index].real,3)
            generation_Q =  round(S_inj[bus_index].imag,3)
            loading = (generation_P**2 + generation_Q**2)**(1/2)*MVA_base/Gen_MVA[k]

        bus_results.append([bus_index+1, bus_label, bus_voltage_mag, bus_voltage_ang, generation_P, generation_Q, loading, load_P, load_Q])

    # Branch flow results
    branch_results = []
    for i in range(len(br_f)):
        from_bus = ind_to_bus[br_f[i] + 1]
        to_bus = ind_to_bus[br_t[i] + 1]
        from_bus_injection_P = V[br_f[i]] * np.conj(Y_from[i].dot(V)).real / MVA_base
        from_bus_injection_Q = V[br_f[i]] * np.conj(Y_from[i].dot(V)).imag / MVA_base
        to_bus_injection_P = V[br_t[i]] * np.conj(Y_to[i].dot(V)).real / MVA_base
        to_bus_injection_Q = V[br_t[i]] * np.conj(Y_to[i].dot(V)).imag / MVA_base

        branch_results.append([i + 1, from_bus, to_bus, from_bus_injection_P, from_bus_injection_Q, to_bus_injection_P, to_bus_injection_Q])

    #show results
    headers_bus = ["Bus Nr", "Label", "Voltage Mag (pu)", "Voltage Ang (deg)", "Generation P (pu)", "Generation Q (pu)","Generation loading", "Load P (pu)", "Load Q (pu)"]
    headers_branch = ["Branch Nr", "From Bus", "To Bus", "From Bus Inject. P (pu)", "From Bus Inject. Q (pu)", "To Bus Inject. P (pu)", "To Bus Inject. Q (pu)"]

    print("=====================================================================")
    print("| Bus results |")
    print("=====================================================================")
    print(tabulate(bus_results, headers=headers_bus, tablefmt="fancy_grid"))
    print("=====================================================================")

    print("\n=========================================================")
    print("| Branch flow |")
    print("=========================================================")
    print(tabulate(branch_results, headers=headers_branch, tablefmt="fancy_grid"))
    print("=========================================================")

    
    return 
