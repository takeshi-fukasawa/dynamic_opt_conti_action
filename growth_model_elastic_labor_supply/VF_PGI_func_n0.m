function [out,other_vars]=VF_PGI_func_n0(V,var2,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec)

   % Written by Takeshi Fukasawa in June 2024, based on the code of Maliar and Maliar (2013)
   
   global lambda_param difference_n0
   
   if Method==-1 %%%%%%%% VF-PGI update
       n0=var2;%initial action
       
       c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
       k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);
       
       vf_coef=X0\V;
       [V_new] =VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                % Recompute value function using the Bellman 
                                % equation

       vf_coef=X0\V; % Coefficients for value function

       EVder=EVder_func(k1,z1,n_nodes,weight_nodes,vf_coef,D);
       difference_n0=FOC_L(EVder,n0,c0,k0,z0,A,alpha,gam,nu,B,beta);

       
       n0_new=n0+lambda_param*difference_n0;

       V_new = kdamp*V_new + (1-kdamp)*V;   
                   % Update V using damping
                      
       out={V_new,n0_new};
        %==================================================================
   end

    other_vars.k1=k1;
    other_vars.c0=c0;
    other_vars.n0=n0;
    other_vars.k0=k0;

return
