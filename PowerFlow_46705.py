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
    
    F= np.concatenate((Delta_P[pv_index],Delta_P[pq_index],Delta_Q[pq_index]),axis=0)

    return F


# 3. the CheckTolerance() function
def CheckTolerance(F,n,err_tol):

    return success

# 4. the generate_Derivatives() function
def generate_Derivatives(Ybus,V):

    return J_ds_dVm,J_dS_dTheta


# 5. the generate_Jacobian() function
def generate_Jacobian(J_dS_dVm,J_dS_dTheta,pv_index,pq_index):

    return J


# 6. the Update_Voltages() function
def Update_Voltages(dx,V,pv_index,pq_index):

    return V



####################################################
#  Displaying the results in the terminal window   #
####################################################
def DisplayResults(V,lnd):

    return
