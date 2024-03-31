def System_violations(V,Ybus,Y_from,Y_to,lnd):
# Inputs:
# V = results from the load flow
# Ybus = the bus admittance matrix used in the load flow
# Y_from,Y_to = tha admittance matrices used to determine the branch flows
# lnd = the LoadNetworkData object for easy access to other model data

#store variables as more convenient names
    br_f=lnd.br_f; br_t=lnd.br_t; #from and to branch indices (FROM LND)
    ind_to_bus=lnd.ind_to_bus; # the ind_to_bus mapping object (FROM LND)
    bus_to_ind=lnd.bus_to_ind; # the bus_to_ind mapping object (FROM LND)
    #br_MVA = lnd.br_MVA # object containing the MVA ratings of the branches (FROM lnd. We have to create this)
    #br_id = lnd.br_id # (you have to update LoadNetworkData for this) (From lnd. We have to create this,but I think is "id_")
    MVA_base = lnd.MVA_base #(I added it)
    bus_labels = lnd.bus_labels #(I added it)
    buscode = lnd.buscode #(I addeed it)
    #MVA_Size=lnd.MVA_Size #(I added it)
    Lines=lnd.Lines #(I added it)
    Gen_MVA=lnd.Gen_MVA #(I added it)
    
#line flows and generators injection....
    S_to = V[br_t]*(Y_to.dot(V)).conj() # the flow in the to end.. 
    S_from = V[br_f]*(Y_from.dot(V)).conj() # the flow in the from end
    S_inj = V*(Ybus.dot(V)).conj() # the injected power
    SLD=lnd.SLD # The defined loads on the PQ busses (Check beacause in loadnetworkdata is defined as SLD)
    S_gen = S_inj + SLD # the generator arrays (I dont know if we use this equation)
    violations = [] #empty list that will store strings describing each violation
    MVAlines = [[Lines[i][-1]] for i in range(len(Lines))] 
    
    
#Check flow in all branches (both ends) and report if limits are violated      
    for i in range(len(Lines)):
        
        S_to = round(V[br_t[i]]*(Y_to.dot(V)).conj()[i],3)
        S_from = round(V[br_f[i]]*(Y_from.dot(V)).conj()[i],3)
        
            
        if abs(S_from)*MVA_base>MVAlines[i][0]: 
            violations.append("Power flow violation in branch %d "%(i+1) + "from bus %d "%(ind_to_bus[br_f[i]]) +  "to bus %d : "%(ind_to_bus[br_t[i]]) + str(round(abs(S_from)/MVAlines[i][0]*MVA_base*100,2)) + "%")
            
        if abs(S_to)*MVA_base>MVAlines[i][0]:
             violations.append("Power flow violation in branch %d "%(i+1) + "from bus %d "%(ind_to_bus[br_t[i]]) +  "to bus %d : "%(ind_to_bus[br_f[i]]) + str(round(abs(S_to)/MVAlines[i][0]*MVA_base*100,2)) + "%")
                  
#Check output of all generators and see if limits are exceeded
    
    for i in range(len(bus_labels)):
        bus_index = i 
        if buscode[bus_index]==2 or buscode[bus_index]==3:
            if abs(S_gen[bus_index])*MVA_base>Gen_MVA[bus_index]:
                violations.append("Power generation violation at bus %d with a value of : " %(ind_to_bus[bus_index])  + str(round(abs(S_gen[bus_index])/Gen_MVA[bus_index]*MVA_base*100,2)) + "%")
    
                
#Check voltages on all busses and see if it remains between 0.9 and 1.1 pu
    N=len(V)
    for i in range (N):
        if abs(V[i])<0.9 :
            violations.append("Voltage violation at bus %d : "%(ind_to_bus[i]) + str(round(abs(V[i]),3)) + " pu < 0.9 pu")
        if abs(V[i])>1.1 :
            violations.append("Voltage violation at bus %d : "%(ind_to_bus[i]) + str(round(abs(V[i]),3)) + " pu < 1.1 pu")
    return violations





