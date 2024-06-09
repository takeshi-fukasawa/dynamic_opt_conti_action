function [out,other_vars]=joint_update_func_L(V,var2,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,spectral_spec)

   % Written by Takeshi Fukasawa in June 2024, based on the code of Maliar and Maliar (2013)
   
   global lambda_param difference_n0
   
   if Method==0 %%%%%%%% VF-PGI update
       n0=var2;%initial action
       
       c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
       k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);
       
       vf_coef=X0\V;
       [V_new] =VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                % Recompute value function using the Bellman 
                                % equation

       vf_coef=X0\V; % Coefficients for value function

       EVder=EVder_func(k1,z1,n_nodes,weight_nodes,vf_coef,D);
       difference_n0=FOC_L_VFI(EVder,n0,c0,k0,z0,A,alpha,gam,nu,B,beta);

       
       n0_new=n0+lambda_param*difference_n0;

       V_new = kdamp*V_new + (1-kdamp)*V;   
                   % Update V using damping
                      
       out={V_new,n0_new};
   

   elseif Method==2
%==================================================================
        % Method 2. Endogenous grid method iterating on value function (EGM-VF)
        %==================================================================            
            k0=var2;
            grid(:,1) = k0;        % Grid points for current capital 
            
            X0 = Polynomial_2d(grid,D);   % Construct polynomial on 
                                          % current state variables


            vf_coef=X0\V; % Coefficients for value function
                       
            k1 = grid_EGM(:,1); % Grid points for next-period capital 
                                % (fixing endogenous grid)
            for j = 1:n_nodes                   
                Xder1_EGM = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1_EGM(:,j) = Xder1_EGM*vf_coef; 
                                % Compute derivative of value function in
                                % the integration nodes    
            end
            
            Wder1 = Vder1_EGM*weight_nodes;   
                                % Compute the expected derivative of next-period 
                                % value function
            for j=1:n_grid      % Solve for labor using eq. (17) in MM (2013)
                n0(j,1)=csolve('Labor_EGM',n0(j,1),[],0.000001,5,nu,gam,alpha,delta,beta,B,A,k1(j,1),z0(j,1),Wder1(j,1));
            end                                                
            
            c0 = (beta*Wder1).^(-1/gam);      
                                % Compute consumption from eq. (5) in MM
                                % (2016) 
            k0 = (B*(1-n0).^-nu./(1-alpha)/A./z0/beta./Wder1).^(1/alpha).*n0;
                                % Compute current capital from eq. (4) in
                                % MM (2013)
                                              
            [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                % Recompute value function using the Bellman 
                                % equation
                                 
            V_new = kdamp*V_new + (1-kdamp)*V;   
                                % Update V using damping

            k0_new=k0;
   
            out={V_new,k0_new};

 
        %==================================================================
   end

    other_vars.k1=k1;
    other_vars.c0=c0;
    other_vars.n0=n0;
    other_vars.k0=k0;

return
