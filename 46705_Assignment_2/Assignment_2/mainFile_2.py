# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
import PowerFlow_46705_2 as pf # importPowerFlowfunctions
import LoadNetworkData_2 as lnd
#import System_violations as SV #agregado
max_iter=30 #Iterationsettings
err_tol=1e-4
#lnd.LoadNetworkData('TestSystem.txt') # LoadtheNetworkdata...
#lnd.LoadNetworkData('Kundur_two_area_system.txt')
lnd.LoadNetworkData('Nordic32_SA.txt')
V,success,n= pf.PowerFlowNewton(lnd.Ybus,lnd.Sbus,lnd.V0,lnd.pv_index,lnd.pq_index,max_iter,err_tol)
# Displayresultsif thepowerflow analysisconverged
if success:
  pf.DisplayResults_and_loading(V,lnd)


#SV.System_violations(V,lnd.Ybus,lnd.Y_fr,lnd.Y_to,lnd) #agregado
