function [out,other_vars]=joint_update_func(V,var2,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,k1_temp,gam,c0,...
   beta,n_nodes,weight_nodes,vf_coef,D,kdamp,n_grid,opts,spectral_spec)

   % Written by Takeshi Fukasawa in March 2024, based on the code of AMMT (2016)
   
   global lambda_param
   
   if Method==8 %% VF-PGI update
       k1=var2;%initial action
       vf_coef=X0\V; % Coefficients for value function
                
       for i = 1:n_grid    % For each grid point 
           difference_i=FOC_VFI(k1(i,1),k0(i,1),z0(i,1),A,alpha,gam,delta,beta,z1(i,:),n_nodes,weight_nodes,vf_coef,D);
           k1_new(i,1)=k1(i,1)+lambda_param*difference_i;
       end
      
       c0 = (1-delta)*k0+A*z0.*k0.^alpha-k1;% Find consumption from budget constraint
       
       [V_new] = VF_Bellman(c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D);
                   % Recompute value function using Bellman equation
    
       V_new = kdamp*V_new + (1-kdamp)*V;   
                   % Update V using damping
                      
       out={V_new,k1_new};
   
       %%%%%%%%%%%%%%%%%%
   elseif Method==3  %% EGM update
            k0=var2;
            grid(:,1) = k0;        % Grid points for current capital 
            
            X0 = Polynomial_2d(grid,D);   % Construct polynomial on 
                                          % current state variables

            vf_coef=X0\V; % Coefficients for value function
            
            k1 = grid_EGM(:,1); % Grid points for next-period capital 
                                % (fixing endogenous grid)
            for j = 1:n_nodes   
                
                X1_EGM    = Polynomial_2d([k1 z1(:,j)],D);                
                V1_EGM(:,j) = X1_EGM*vf_coef; % Compute value function in 
                                              % the integration nodes    

                Xder1_EGM = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1_EGM(:,j) = Xder1_EGM*vf_coef; % Compute derivative 
                                                    % of value function 
                                                    % in the integration nodes    
            end
            c0 = (beta*Vder1_EGM*weight_nodes).^(-1/gam);
            
            for i = 1:n_grid       % For each grid point 
                k0(i) = fsolve('EGM_BC',k0(i),opts,k1(i),z0(i),c0(i),A,alpha,delta);
                                   % Solve for current capital that satisfies 
                                   % budget constraint given next-period capital
            end
            

            if gam==1
                V_new = log(c0)+beta*V1_EGM*weight_nodes; 
                                   % Bellman equation if gam=1        
            else
                V_new = (c0.^(1-gam)-1)/(1-gam)+beta*V1_EGM*weight_nodes;
            end                    % Bellman equation otherwise
                     
            V_new = kdamp*V_new + (1-kdamp)*V;   
                                   % Update V using damping
            k0_new=k0;
   
            out={V_new,k0_new};
   
        %==================================================================
   end

    other_vars.k1=k1;
    other_vars.c0=c0;
    other_vars.k0=k0;
    other_vars.grid=grid;
    other_vars.X0=X0;
    other_vars.V=V;

return

