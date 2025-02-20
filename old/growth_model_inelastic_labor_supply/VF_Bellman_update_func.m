function [out,other_vars]=VF_Bellman_update_func(V,...
    X0,c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D,kdamp)

        vf_coef = X0\V; 
        [V_new] = VF_Bellman(c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D);
        V_new = kdamp*V_new + (1-kdamp)*V; 
 
      out={V_new};
      other_vars=[];
return
