function [out,other_vars]=V_update_func(V,X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D,kdamp)


vf_coef=X0\V;
            
    [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                 % Recompute value function using 
                                 % the Bellman equation

    V_new=kdamp*V_new+(1-kdamp)*V; % Update V using damping

out={V_new};
other_vars=[];

end
