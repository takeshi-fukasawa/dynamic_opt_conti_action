function [out,other_vars]=VF_Bellman_L_update_func(V_init,...
    X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef_temp,D,kdamp)

        vf_coef = X0\V_init; 
        [V_new] =VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
        V_new = kdamp*V_new + (1-kdamp)*V_init; 
 
      out={V_new};
      other_vars=[];
return