def apply_contingency_to_Y_matrices(Ybus,Yfr,Yto,fr_ind,to_ind,br_ind,Ybr_mat):
    # input:
    # The original admitance matirces: Ybus,Yfr,Yto
    # The from and to end indices for the branch (fr_ind, to_ind)
    # The indice for where the branch is in the branch list (br_ind)
    # The 2x2 admittance matrix for the branch Ybr_mat
    ##########################################################
    # This is important, you must copy the original matrices
    Ybus_mod = Ybus.copy() # otherwise you will change the original Ybus matrix
    Yfr_mod = Yfr.copy() # when ever you make changes to Ybus_mod
    Yto_mod = Yto.copy() # using the .copy() function avoids this
    ##################################################################################
    #
    # YOUR CODE COMES HERE:
    #
    # 1. Remove the branch from the Ybus_mod matrix

    Ybus_mod[fr_ind,fr_ind] -= Ybr_mat[0,0]
    Ybus_mod[fr_ind,to_ind] -= Ybr_mat[0,1]
    Ybus_mod[to_ind,fr_ind] -= Ybr_mat[1,0]
    Ybus_mod[to_ind,to_ind] -= Ybr_mat[1,1]
    
    # 2. Remove the branch from the Yto and Yfr matrices
    
    Yfr_mod[br_ind,fr_ind] -= Ybr_mat[0,0]       
    Yfr_mod[br_ind,to_ind] -= Ybr_mat[0,1]
    Yto_mod[br_ind,to_ind] -= Ybr_mat[1,1]       
    Yto_mod[br_ind,fr_ind] -= Ybr_mat[1,0]
    
    
    ####################################################################################


    
    
    
    
    return Ybus_mod,Yfr_mod,Yto_mod


#breakpoint()