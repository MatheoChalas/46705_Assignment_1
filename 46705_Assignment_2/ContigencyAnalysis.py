# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:04:56 2024

@author: mathe
"""

import numpy as np
import PowerFlow_46705 as pf #importPowerFlow functions
import LoadNetworkData_V2023 as lnd # loadthenetworkdatato globalvariables
max_iter = 30 #Iterationsettings
err_tol = 1*10**(-4)
# LoadtheNetworkdata...
filename = "./TestSystem4SA.txt"
lnd.LoadNetworkData(filename) # makesYbusavailableas lnd.Ybusetc.
#%%
#######################################################################
# Part I: Studythebase caseanddisplayresults(with% loading)#
#######################################################################
V,success,n = pf.PowerFlowNewton(lnd.Ybus,lnd.Sbus,lnd.V0,lnd.pv_index,
lnd.pq_index,max_iter,err_tol)
if success: # Displayresultsif thepowerflow analysisconverged
    pf.DisplayResults_and_loading(V,lnd)

#%%
######################################################################
# Part II: Simplifiedcontingencyanalysis(onlybranchoutages) #
######################################################################
print('*'*50)
print('* ContingencyAnalysis *')
print('*'*50)

for i in range(len(lnd.br_f)): #sweepoverbranches
    fr_ind = lnd.br_f[i]
    to_ind = lnd.br_t[i]
    br_ind = i
    Ybr_mat = lnd.br_Ymat[i]
    Ybus_mod,Yfr_mod,Yto_mod = \
        apply_contingency_to_Y_matrices(lnd.Ybus,lnd.Y_fr,lnd.Y_to,fr_ind,to_ind,br_ind,Ybr_mat)

    str_status = '-'*63 + '\nTrippingof branch{:}(bus{:}- bus{:})'.format(i+1,lnd.ind_to_bus[fr_ind] ,lnd.ind_to_bus[to_ind])

    try: #try theloadflow,if it fails,displaymessage
        V,success,n = pf.PowerFlowNewton(Ybus_mod,lnd.Sbus,lnd.V0,lnd.pv_index,lnd.pq_index,max_iter,err_tol,print_progress=False)
    except:
        str_ = '--> LoadFlowerror(Jacobian) whenbranch{:}(bus {:}- bus{:})is tripped'
        str_ = str_.format(i+1,lnd.ind_to_bus[fr_ind] ,lnd.ind_to_bus[to_ind])
        print(str_status+ ' [CONVERGENCEISSUES!]')
        print(str_)
    else:
        if success: # Displayresultsif thepowerflow analysisconverged
            violations = System_violations(V,Ybus_mod,Yfr_mod,Yto_mod,lnd)
            if not violations: # no violations, printstatusandmoveon
                print(str_status + '[OK!]')
            else: # ifviolation, displaythem
                print(str_status + '[Violations!]')
                for str_ in violations:
                    print(str_)
        else: #noconvergence...
            str_ = '-->No load-flowconvergencewhen branch{:}(bus{:}-bus {:})is tripped'
            str_ = str_.format(i+1,lnd.ind_to_bus[fr_ind] ,lnd.ind_to_bus[to_ind])
            print(str_status + ' [CONVERGENCEISSUES!]')
            print(str_)