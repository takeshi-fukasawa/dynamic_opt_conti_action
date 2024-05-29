function [out,other_vars]=update_func(input_cell,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,k1,gam,c0,...
   beta,n_nodes,weight_nodes,vf_coef,D,kdamp,n_grid,opts,spectral_spec)
  
  % Based on AMMT (2016)
  % Modified by Takeshi Fukasawa in March 2024

  global V X0_temp

if Method<=3 & Method>=1
    V=input_cell;
elseif Method==4 | Method==5 | Method==7
    k1=input_cell;
elseif Method==6
   Vder0=input_cell;
end


%%%%%%%% Based on the original code of Main_7_Methods.m %%%%

        %  Method 1. Envelope condition method iterating on value function (ECM-VF)
        %==================================================================
        if Method==1        
            vf_coef=X0\V; % Coefficients for value function

            Vder0 = X0der*vf_coef;      % Compute derivative of value function
            u0 = Vder0./(1-delta+A*alpha*z0.*k0.^(alpha-1));
                                        % Computer marginal utility from
                                        % envelope condition
            c0 = u0.^(-1/gam);          % Find consumption   
            k1 = (1-delta)*k0+A*z0.*k0.^alpha-c0;  
                                        % Find capital from budget constraint
            [V_new] = VF_Bellman(c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D);
                                        % Recompute value function using 
                                        % Bellman equation
            V_new=kdamp*V_new+(1-kdamp)*V; % Update V using damping

        %==================================================================

        
        % Method 2. Conventional value function iteration (VFI)
        %==================================================================
        elseif Method==2;            
            vf_coef=X0\V; % Coefficients for value function
            
            for i = 1:n_grid       % For each grid point 
                k1(i) = fsolve('FOC_VFI',k1(i),opts,k0(i),z0(i),A,alpha,gam,delta,beta,z1(i,:),n_nodes,weight_nodes,vf_coef,D);
                                   % Solve for capital that satisfies FOC
            end
            c0 = (1-delta)*k0+A*z0.*k0.^alpha-k1; 
                                   % Find consumption from budget constraint
            [V_new] = VF_Bellman(c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D);
                                   % Recompute value function using Bellman 
                                   % equation

            V_new = kdamp*V_new + (1-kdamp)*V;   
                                   % Update V using damping
        %==================================================================
                                                             
                                   
        % Method 3. Endogenous grid method (EGM)
        %==================================================================                                   
        elseif Method==3;         
            X0=X0_temp;
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
                     
            grid(:,1) = k0;        % Grid points for current capital 
            
            X0 = Polynomial_2d(grid,D);   % Construct polynomial on 
                                          % current state variables
            X0_temp=X0;%%%%

            V_new = kdamp*V_new + (1-kdamp)*V;   
                                   % Update V using damping
        %==================================================================
                                                            
                                                           
        % Method 4. Policy function iteration via envelope condition (ECM-PI)
        %==================================================================                                                                     
        elseif Method==4 
       
                                 % given policy function
            c0 = (1-delta)*k0+A*z0.*k0.^alpha-k1;    
                                 % Compute consumption
           
            % Solve for value function from Bellman equation
            spec.TOL=1e-6;
            spec.norm_spec=10;
            spec.update_spec=spectral_spec; 
            [output_spectral,other_vars,iter_info]=...
                spectral_func(@VF_Bellman_update_func,spec,{V},...
                   X0,c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

           V=output_spectral{1};
           vf_coef = X0\V;     % Coefficients for value function
                         
            Vder0 = X0der*vf_coef;      % Compute derivative of value function
            u0 = Vder0./(1-delta+A*alpha*z0.*k0.^(alpha-1));
                                        % Compute marginal utility from
                                        % envelope condition
            c0 = u0.^(-1/gam);          % Find consumption   
            k1_new = (1-delta)*k0+A*z0.*k0.^alpha-c0;  
                                        % Find capital from budget constraint
            k1_new = kdamp*k1_new + (1-kdamp)*k1;               
                                        % Update k1 using 
                                        % damping
        %==================================================================
  
        
        
        % Method 5. Conventional policy function iteration via FOC (PI)
        %==================================================================                                                               
        elseif Method==5  
                          
            c0 = (1-delta)*k0+A*z0.*k0.^alpha-k1;    
                                 % Compute consumption

            spec.TOL=1e-6;
            spec.norm_spec=10; 
            [output_spectral,other_vars,iter_info]=...
                spectral_func(@VF_Bellman_update_func,spec,{V},...
                   X0,c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

           V=output_spectral{1};
           vf_coef = X0\V;     % Coefficients for value function
          
            % Recompute the capital policy function using FOC
            for i = 1:n_grid     % For each grid point 
                k1_new(i,1) = fsolve('FOC_VFI',k1(i),opts,k0(i),z0(i),A,alpha,gam,delta,beta,z1(i,:),n_nodes,weight_nodes,vf_coef,D);
                                 % Solve for capital that satisfies FOC
            end
           k1_new = kdamp*k1_new + (1-kdamp)*k1;               
                                        % Update k1 using 
                                        % damping
        %==================================================================
        
        
        
        % Method 6. Envelope condition method iterating on derivative of 
        % value function (ECM-DVF)
        %==================================================================
        elseif Method==6  
            vf_coef=X0der\Vder0; % Coefficients for value function
        
            u0 = Vder0./(1-delta+A*alpha*z0.*k0.^(alpha-1));
                                        % Compute marginal utility from
                                        % envelope condition
            c0 = u0.^(-1/gam);          % Find consumption   
            k1 = (1-delta)*k0+A*z0.*k0.^alpha-c0;  
                                        % Find capital from budget constraint
                                        
            for j = 1:n_nodes           
                X1der = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1(:,j) = X1der*vf_coef; % Compute derivative of value
                                            % function in the integration nodes
            end
            
            Vder0_new = (1-delta+A*alpha*z0.*k0.^(alpha-1)).*(beta*Vder1*weight_nodes);
                                        % Recompute derivative of value 
                                        % function in the grid points
            warning('off')              % Some polynomial terms are zero for
                                        % the derivative and system is 
                                        % underdetermined. Least square 
                                        % problem is still correctly 
                                        % processed by a truncated QR method
                                        % but the system produces a warning
            Vder0_new = kdamp*Vder0_new + (1-kdamp)*Vder0;   
                                   % Update Vder0 using damping
        %==================================================================
        
        
        
        % Method 7. Conventional Euler equation method (EE)
        %==================================================================               
        elseif Method==7  
            K_coef = X0\k1;     % New coefficients for policy function 
            
            for j = 1:n_nodes    
                X1 = Polynomial_2d([k1 z1(:,j)],D);            
                k2(:,j) = X1*K_coef; % Compute capital in the integration nodes
            end
            
            k1_dupl = k1*ones(1,n_nodes);
                     % Duplicate k1 n_nodes times to create a matrix with 
                     % n_nodes identical columns;   
            c1  = (1-delta)*k1_dupl+A*z1.*k1_dupl.^alpha-k2;
                     % Compute consumption in all integration nodes
            c0  =  (beta*(c1.^(-gam).*(1-delta+alpha*A*k1_dupl.^(alpha-1).*z1))*weight_nodes).^(-1/gam);
                     % Compute current consumption using Euler equation                    
            k1_new =(1-delta)*k0+A*z0.*k0.^alpha-c0;    
                     % Compute new capital on the grid
            k1_new = kdamp*k1_new + (1-kdamp)*k1;               
                    % Update the coefficients using damping      
        %==================================================================
                          
        end               
        
%%%%%%%%
if Method<=3 & Method>=1
    out={V_new};
elseif Method==4 | Method==5 | Method==7
    out={k1_new};
elseif Method==6
    out={Vder0_new};
end

other_vars.k1=k1;
other_vars.c0=c0;
other_vars.k0=k0;
other_vars.grid=grid;
other_vars.X0=X0;
other_vars.V=V;

return

