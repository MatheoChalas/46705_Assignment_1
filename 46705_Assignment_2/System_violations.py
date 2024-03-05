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
    br_MVA = lnd.br_MVA # object containing the MVA ratings of the branches (FROM lnd. We have to create this)
    br_id = lnd.br_id # (you have to update LoadNetworkData for this) (From lnd. We have to create this,but I think is "id_")
    MVA_base = lnd.MVA_base #(I added it)
    bus_labels = lnd.bus_labels #(I added it)
    buscode = lnd.buscode #(I addeed it)
    MVA_Size=lnd.MVA_Size #(I added it)
    
#line flows and generators injection....
    S_to = V[br_t]*(Y_to.dot(V)).conj() # the flow in the to end.. 
    S_from = V[br_f]*(Y_from.dot(V)).conj() # the flow in the from end
    S_inj = V*(Ybus.dot(V)).conj() # the injected power
    SLD=lnd.S_LD # The defined loads on the PQ busses (Check beacause in loadnetworkdata is defined as SLD)
    S_gen = S_inj + SLD # the generator arrays (I dont know if we use this equation)
    violations = [] #empty list that will store strings describing each violation
    
    
#Check flow in all branches (both ends) and report if limits are violated      
    for i in range(len(br_f)):
        
        S_to = round(V[br_t[i]]*(Y_to.dot(V)).conj()[i],3)
        S_from = round(V[br_f[i]]*(Y_from.dot(V)).conj()[i],3)
        
        if abs(S_to)*MVA_base>br_MVA[i]:
            violations.append(S_to)
            
        if abs(S_from)*MVA_base>br_MVA[i]:
            violations.append(S_from)
               
#Check output of all generators and see if limits are exceeded
    for i in range(len(bus_labels)):
        bus_index = i
        if buscode[bus_index]==2:
            if abs(S_inj[bus_index])*MVA_base>MVA_Size[bus_index]:
                violations.append(S_inj[bus_index])
                
#Check voltages on all busses and see if it remains between 0.9 and 1.1 pu
    N=len(V)
    for i in range (N):
        if abs(V[i])<0.9 or abs(V[i])>1.1:
            violations.append(V[i])
            
    return violations





