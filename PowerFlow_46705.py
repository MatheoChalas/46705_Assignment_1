"""
46705 - Power Grid Analysis
This file contains the definitions of the functions needed to
carry out Power Flow calculations in python.

How to carry out Power Flow in a new *.py file? 
See the example in table 1 in the assignment text
"""

import numpy as np


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
     while (notsuccess) and (n < max_iter):
         n += 1 # Updatecounter
         # ComputederivativesandgeneratetheJacobianmatrix
         J_dS_dVm,J_dS_dTheta = generate_Derivatives(Ybus,V)
         J = generate_Jacobian(J_dS_dVm,J_dS_dTheta,pv_index,pq_index)
         # Computetheupdatestep
         dx = np.linalg.solve(J,F)
         # Updatevoltagesand checkiftoleranceisnow reached
         V =Update_Voltages(dx,V,pv_index,pq_index)
         F = calculate_F(Ybus,Sbus,V,pv_index,pq_index)
         success = CheckTolerance(F,n,err_tol)
     
     if success: #printoutmessageconcerningwetherthe powerflowconvergedornot
         print('TheNewtonRapsonPowerFlow Convergedin %diterations!'% (n,))
     else:
         print('No Convergence!!!\n Stoppedafter%diterationswithoutsolution...'% (n,))
     return V,success,n


# 2. the calculate_F() function
def calculate_F(Ybus,Sbus,V,pv_index,pq_index):
    Delta_S= Sbus- V * (Ybus.dot(V)).conj()

    Delta_P[pv_index]=Delta_S[pv_index].real
    Delta_P[pq_index]=Delta_S[pq_index].real
    Delta_q[pq_index]=Delta_S[pq_index].imag
    
    F= np.concatenate((Delta_P[pv_index],Delta_P[pq_index],Delta_Q[pq_index]),axis=0)

    return F


# 3. the CheckTolerance() function
def CheckTolerance(F,n,err_tol):
    normF = np.linalg.norm(F,np.inf)
    if normF<err_tol:
        success=1
        print(normF,n)
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
def DisplayResults(V,lnd):

    return
