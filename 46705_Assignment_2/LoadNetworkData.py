# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:00:27 2024

@author: mathe
"""

 #################################################
 # New code to be included in LoadNetworkData() #
 #################################################
global Lines, Trans, Gen_MVA #new variables as global variables.
Lines = []
Trans = []
Gen_MVA = np.zeros(N) #keep track of generators MVA size (bus indices used)
#########################################################
# modifications when dealing with transmission lines... #
#########################################################
for line in line_data:
    bus_fr, bus_to, id_, R, X, B, X2, X0, MVA_rate = line #unpack the values
    Lines.append([bus_fr, bus_to, id_,Y_se,Y_sh_2,MVA_rate]

                 ###############################################
# ..... the remaining code comes here
#(for the Ybus matrices ....)
##############################################
#########################################################
# modifications when dealing with transformers...
#
#########################################################

for line in tran_data:
    bus_fr, bus_to, id_, R,X,n,ang1,fr_co, to_co, X2, X0, MVA_rate = line #unpack values
    Trans.append([bus_fr, bus_to, id_, Yps_mat, MVA_rate ])

###############################################
# ..... the remaining code comes here
#(for the Ybus matrices ....)
#############################################