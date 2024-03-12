# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
import PowerFlow_46705_2 as pf # importPowerFlowfunctions
import LoadNetworkData_2 as lnd
import System_violations as sv #agregado
import apply_contingency_to_Y_matrices as ac

max_iter=30 #Iterationsettings
err_tol=1e-4
#lnd.LoadNetworkData('TestSystem4SA.txt')
lnd.LoadNetworkData('Nordic32_SA.txt')
V,success,n= pf.PowerFlowNewton(lnd.Ybus,lnd.Sbus,lnd.V0,lnd.pv_index,lnd.pq_index,max_iter,err_tol)
# Displayresultsif thepowerflow analysisconverged
if success:
  pf.DisplayResults_and_loading(V,lnd)


print('*'*50 )
print('*    Contingency Analysis    *')
print('*'*50 )

for i in range(len(lnd.br_f) ) : #sweep over branches
    fr_ind = lnd.br_f[i]
    to_ind = lnd.br_t[i]
    br_ind = i
    Ybr_mat = lnd.br_Ymat[i]
    Ybus_mod , Yfr_mod , Yto_mod =ac.apply_contingency_to_Y_matrices( lnd.Ybus ,lnd.Y_fr,lnd.Y_to,fr_ind,to_ind,br_ind,Ybr_mat)

    str_status = '-'*63 + '\nTripping of branch {:} (bus {:} - bus{:})'. format(i+ 1,lnd.ind_to_bus[fr_ind] ,lnd.ind_to_bus[to_ind])

    try : #try the load flow, if it fails, display message
        V , success , n = pf.PowerFlowNewton(Ybus_mod , lnd . Sbus , lnd . V0 , lnd . pv_index , lnd . pq_index ,max_iter , err_tol , print_progress=False)
    except :
        str_ = ' --> Load Flow error (Jacobian) when branch {:} (bus {:} - bus {:}) is tripped'
        str_ = str_ . format(i+1 ,lnd . ind_to_bus[fr_ind] , lnd . ind_to_bus[to_ind ] )
        print(str_status+ ' [CONVERGENCE ISSUES!]')
        print(str_)
    else :
        if success : # Display results if the power flow analysis converged
            violations =sv.System_violations( V , Ybus_mod , Yfr_mod , Yto_mod , lnd)
            if not violations : # no violations, print status and move on
                print(str_status + ' [OK!]')
            else : # if violation, display them
                print(str_status + ' [Violations!]')
                for str_ in violations:
                    print(str_)
        else : #no convergence...
            str_ = ' --> No load-flow convergence when branch {:} (bus {:} - bus {:}) is tripped' 
            str_ = str_ . format(i+1 ,lnd . ind_to_bus[fr_ind ] , lnd . ind_to_bus[to_ind ] )
            print(str_status + ' [CONVERGENCE ISSUES!]')
            print(str_) 

