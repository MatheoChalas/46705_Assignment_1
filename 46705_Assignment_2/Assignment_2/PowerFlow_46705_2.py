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
    
    #Computation of the apparent power
    Delta_S= Sbus- V * (Ybus.dot(V)).conj()
    
    #Computation of the active power and reactive power
    Delta_P=Delta_S.real
    Delta_Q=Delta_S.imag
    
    #Computation of the matrix F which is the mismatch vector corresponding to the mismatch between the specified values and the calculated ones
    F= np.concatenate((Delta_P[pv_index],Delta_P[pq_index],Delta_Q[pq_index]),axis=0)

    return F


# 3. the CheckTolerance() function
def CheckTolerance(F,n,err_tol):
    #Computation of the maximal mismatch
    normF = np.linalg.norm(F,np.inf)
    
    #Display of the first line of the Test Tolerance
    if n==0:
        print("Check Tolerance :")
    print("Absolute value of the greatest mismatch : %f" %normF,"Iteration number : %d" %n)
    
    #Display of the last line of the Test Tolerance and giving to success the value of 1 if the mismatch is lower than the error tolerated
    if normF<err_tol:
        success=1
        print("End of Check Tolerance")
        print("------------------------------------")
    
    #If the maximal mismatch is not below the error tolerated it gives the value of 0 to success 
    else :
        success=0
        
    #Return the value of success to check if we reached an acceptable mismatch for F
    return success

# 4. the generate_Derivatives() function
def generate_Derivatives(Ybus,V):
    #Computation of the derivatives with respect to the voltage magnitude
    J_ds_dVm=np.diag(V/np.absolute(V)).dot(np.diag((Ybus.dot(V)).conj()))+ np.diag(V).dot(Ybus.dot(np.diag(V/np.absolute(V))).conj())
    
    #Computation of the derivatives with respect to the voltage angles
    J_dS_dTheta = 1j*np.diag(V).dot((np.diag(Ybus.dot(V))-Ybus.dot(np.diag(V))).conj())

    #Return the derivatives
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
    
    N1 = 0; N2 = len(pv_index) # dx[N1:N2]-ang.onthe pvbuses
    N3 = N2; N4 = N3 + len(pq_index)# dx[N3:N4]-ang.onthe pqbuses
    N5 = N4; N6 = N5 + len(pq_index)# dx[N5:N6]-mag.onthe pqbuses
    
    #Update of the voltages
    Theta =np.angle(V);Vm =np.absolute(V)
    #Test if the system has pv busses 
    if len(pv_index)>0:
        #Summation of the dx corresponding to the pv busses with the voltage angles
        Theta[pv_index]+= dx[N1:N2]
    if len(pq_index)>0:
        #Summation of the dx corresponding to the pq busses with the voltage angles first and then the voltage magnitudes
        Theta[pq_index]+= dx[N3:N4]
        Vm[pq_index]+=dx[N5:N6]
    
    #Computation of the Voltages
    V = Vm * np.exp(1j*Theta)
    
    return V



####################################################
#  Displaying the results in the terminal window   #
####################################################
from tabulate import tabulate
def DisplayResults_and_loading(V,lnd):
    
    #Importation of all the values from the file lnd that we need to display the results
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
    Gen_MVA=lnd.Gen_MVA
    Lines=lnd.Lines
    
    
    #Computation of the injected apparent powers
    S_inj = V*(Ybus.dot(V)).conj()
    bus_results = []
    
    
    pq_index = np.where(buscode == 1)[0] # Find indices for all PQ-busses
    pv_index = np.where(buscode == 2)[0] # Find indices for all PV-busses
    ref = np.where(buscode == 3)[0] # Find index for ref bus
    
    #Extraction of the needed values for all the busses 
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
        loading="-"
        loading_from="-"
        loading_to="-"
        
        
        #Differenciation of the loads and generators
        if buscode[bus_index]==1:
            load_P = round(SLD[bus_index].real,3)
            load_Q = round(SLD[bus_index].imag,3)
            
        else:
            generation_P = round(S_inj[bus_index].real+SLD[bus_index].real,3)
            generation_Q =  round(S_inj[bus_index].imag+SLD[bus_index].imag,3)
            loading="%.2f" %round(((((generation_P**2)+(generation_Q**2))**0.5)/Gen_MVA[bus_index])*MVA_base*100,2) +"%"
            load_P = SLD[bus_index].real
            load_Q = SLD[bus_index].imag
            

        bus_results.append([bus_index+1, bus_label, bus_voltage_mag, bus_voltage_ang, generation_P, generation_Q,loading,load_P, load_Q])

    # Branch flow results
    branch_results = []
    MVAlines = [[Lines[i][-1]] for i in range(len(Lines))] #Extraction of the Power Ratings from matrix Lines
    
    #Extraction of the needed values for all the busses 
    for i in range(len(Lines)):
       
        from_bus = ind_to_bus[br_f[i]]
        to_bus = ind_to_bus[br_t[i]]
        
        #Computation of the apparent powers flowing in both directions
        S_to = round(V[br_t[i]]*(Y_to.dot(V)).conj()[i],3)
        S_from = round(V[br_f[i]]*(Y_from.dot(V)).conj()[i],3)
        
        from_bus_injection_P = round(S_from.real ,3)
        from_bus_injection_Q = round(S_from.imag ,3)
        to_bus_injection_P = round(S_to.real,3)
        to_bus_injection_Q = round(S_to.imag,3)
        loading_from="%.2f" %round(((((from_bus_injection_P**2)+(from_bus_injection_Q**2))**0.5)/MVAlines[i][0])*MVA_base*100,2) +"%"
        loading_to="%.2f" %round(((((to_bus_injection_P**2)+(to_bus_injection_Q**2))**0.5)/MVAlines[i][0])*MVA_base*100,2) +"%"
        
        
        #Add the wanted data to the branch results list
        branch_results.append([i + 1, from_bus, to_bus, from_bus_injection_P, from_bus_injection_Q,loading_from, to_bus_injection_P, to_bus_injection_Q,loading_to])
        i+=1
    #show results
    
    #Define the headers of the results
    headers_bus = ["Bus Nr", "Label", "Voltage Mag (pu)", "Voltage Ang (deg)", "Generation P (pu)", "Generation Q (pu)","Loading","Load P (pu)", "Load Q (pu)"]
    headers_branch = ["Branch Nr", "From Bus", "To Bus", "From Bus Inject. P (pu)", "From Bus Inject. Q (pu)","From Bus Loading", "To Bus Inject. P (pu)", "To Bus Inject. Q (pu)","To Bus Loading"]

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
