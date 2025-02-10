function [out,other_vars]=update_func_L(input_cell,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec)
  
  % Based on Maliar and Maliar (2013)
  % Modified by Takeshi Fukasawa in June 2024

   global geval_total feval_V_total n0
   global OPI_param
   global PI_linear_eq_sol_method_spec ECM_spec analytical_EE_spec


if Method==1 | Method==0 | Method==4
    V=input_cell;
elseif Method==2 | Method==8 | Method==9
    vf_coef=input_cell;
end


if Method==0
%==================================================================
% Method 0. Value function iteration (VFI)
%==================================================================
     vf_coef=X0\V;

      %%% Use global n0 as initial values => hot start
     for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
          [n0(j,1),exitflag,n_iter,geval]=csolve('FOC_L_VFI',n0(j,1),[],0.000001,5,...
              k0(j,1),z0(j,1),A,alpha,gam,nu,B,beta,delta,...
              z1(j,:),n_nodes,weight_nodes,vf_coef,D);
          geval_total=geval_total+geval;   
     end
                            
     c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
     k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);% Compute next-period capital using budget 
                     % constraint (2) in MM(2013)

    [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                 % Recompute value function using 
                                 % the Bellman equation


    V_new=kdamp*V_new+(1-kdamp)*V; % Update V using damping
                           

elseif Method==1        
%==================================================================
% Method 1. Envelope condition method iterating on value function (ECM-VF)
%==================================================================
     vf_coef=X0\V;
     Vder0 = X0der*vf_coef;   % Compute the derivative of value function
       
      %%% Use global n0 as initial values => hot start       
     for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
          [n0(j,1),exitflag,n_iter,geval]=csolve('Labor_ECM',n0(j,1),[],0.000001,5,nu,alpha,delta,B,A,k0(j,1),z0(j,1),Vder0(j,1));  
          geval_total=geval_total+geval; 
     end
                            
     c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
     k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);% Compute next-period capital using budget 
                     % constraint (2) in MM(2013)

    [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                 % Recompute value function using 
                                 % the Bellman equation

    V_new=kdamp*V_new+(1-kdamp)*V; % Update V using damping
  
    elseif Method==2
        %==================================================================
        % Method 2. Endogenous grid method iterating on value function (EGM-VF)
        %==================================================================            
        
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

      %%% Use global n0 as initial values => hot start
        for j=1:n_grid      % Solve for labor using eq. (17) in MM (2013)
            [n0(j,1),exitflag,n_iter,geval]=csolve('Labor_EGM',n0(j,1),[],0.000001,5,nu,gam,alpha,delta,beta,B,A,k1(j,1),z0(j,1),Wder1(j,1));
          geval_total=geval_total+geval;
        end                                                

        c0 = (beta*Wder1).^(-1/gam);      
                % Compute consumption from eq. (5) in MM
                % (2016) 
        k0 = (B*(1-n0).^-nu./(1-alpha)/A./z0/beta./Wder1).^(1/alpha).*n0;
                % Compute current capital from eq. (4) in
                % MM (2013)
         
        grid(:,1) = k0;        % Grid points for current capital 
            
        X0 = Polynomial_2d(grid,D);   % Construct polynomial on 
                                          % current state variables

        if OPI_param==1          
            [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                % Recompute value function using the Bellman 
                % equation
           feval_V_total=feval_V_total+n_grid;

        else%OPI_param>=2

            V=X0*vf_coef;
            [V_new,feval_V_total]=policy_eval_func(V,n0,c0,k1,...
            X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,gam,...
            nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec);

        end

         vf_coef_new=X0\V_new;

         vf_coef_new = kdamp*vf_coef_new + (1-kdamp)*vf_coef;   
                % Update vf_coef using damping

elseif Method==3
%==================================================================
% Method 3. Euler equation method
%==================================================================

    n0=input_cell;

    c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
    k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);
    
    n_coef=X0\n0;

    dEV_dk1=0;
    for j = 1:n_nodes    % Sum up the expected derivative of value
            % function across the integration nodes
            X1_j = Polynomial_2d([k1 z1(:,j)],D);  
            % Construct the derivatives of basis functions 
            % of polynomial approximating value function at
            % integration node j

            n1_j=X1_j*n_coef;

            numer=B*(1-n1_j).^(-nu);
            denom=(1-alpha).*A.*z1(:,j).*(k1.^alpha).*(n1_j.^(-alpha));

            du1_dc1=numer./denom;

            dEV_dk1=dEV_dk1+du1_dc1.*(1-delta+z1(:,j).*A.*alpha.*(k1.^(alpha-1)).*(n1_j.^(1-alpha))).*weight_nodes(j);
    end

    dk1_dn0=z0.*A.*(k0.^alpha).*(1-alpha).*(n0.^(-alpha));
    dEV_dn0=dk1_dn0.*dEV_dk1;

    if analytical_EE_spec==1% Analytical sol
        numer=B;
        denom=dEV_dn0*beta;
        n0_new=1-(numer./denom).^(1/nu);
    else % Numerically solve sol
        for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
          [n0_new(j,1),exitflag,n_iter,geval]=csolve('EE_resid_func',...
              n0(j,1),[],0.000001,5,...
              dEV_dn0(j),B,beta,nu);  
        end
    end

       out={n0_new};

elseif Method==4
     
%==================================================================
% Method 4. Policy iteration (PI) directly updating V
%==================================================================

     vf_coef=X0\V;
     
     %%% Use global n0 as initial values => hot start  

     if ECM_spec==0
         for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
              [n0(j,1),exitflag,n_iter,geval]=csolve('FOC_L_VFI',n0(j,1),[],0.000001,5,...
                  k0(j,1),z0(j,1),A,alpha,gam,nu,B,beta,delta,...
                  z1(j,:),n_nodes,weight_nodes,vf_coef,D);
              geval_total=geval_total+geval;   
         end
     
     else%%ECM_spec==1
         Vder0 = X0der*vf_coef;   % Compute the derivative of value function
           
          %%% Use global n0 as initial values => hot start       
         for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
              [n0(j,1),exitflag,n_iter,geval]=csolve('Labor_ECM',n0(j,1),[],0.000001,5,nu,alpha,delta,B,A,k0(j,1),z0(j,1),Vder0(j,1));  
              geval_total=geval_total+geval; 
         end
     end

     c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
     k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);% Compute next-period capital using budget 
                     % constraint (2) in MM(2013)

    %% Solve for value function
   [V_new,feval_V_total]=policy_eval_func(V,n0,c0,k1,...
   X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec);

       
   elseif Method==8 
       %==================================================================
        % Method 8. Envelope condition method iterating on derivative of value 
        %           function (ECM-DVF)
        %==================================================================
              
            Vder0=X0der*vf_coef;
        
            for j=1:n_grid           % Solve for labor using eq. (18) in MM (2013)
                [n0(j,1),exitflag,n_iter,geval]=csolve('Labor_ECM',n0(j,1),[],0.000001,5,nu,alpha,delta,B,A,k0(j,1),z0(j,1),Vder0(j,1));
                geval_total=geval_total+geval;      
            end
   
         
           c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
           k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);% Compute next-period capital using budget 
                     % constraint (2) in MM(2013)

            for j = 1:n_nodes           
                X1der = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1(:,j) = X1der*vf_coef; 
                                     % Compute the derivative of value
                                     % function in the integration nodes
            end

            Wder1 = Vder1*weight_nodes;   
                                     % Compute expected derivative  of
                                     % next-period value function
         
            Vder0_new = (1-delta+A*alpha*z0.*k0.^(alpha-1).*n0.^(1-alpha)).*(beta*Wder1);
                                     % Recompute the derivative of value 
                                     % function in the grid points
            
            feval_V_total=feval_V_total+n_grid;

            vf_coef_new = X0der\Vder0_new;

        elseif Method==9;         
        %==================================================================
        %  Method 9. Endogenous grid method iterating on derivative of value 
        %            function (EGM-DVF)
        %==================================================================
            k1 = grid_EGM(:,1);      % Grid points for next-period capital 
                                     % (fixing endogenous grid)
            for j = 1:n_nodes              
                Xder1_EGM = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1_EGM(:,j) = Xder1_EGM*vf_coef; 
                                     % Compute the derivative of value function 
                                     % in the integration nodes    
            end
            
            Wder1 = Vder1_EGM*weight_nodes;   
                                     % Compute expected derivative  of
                                     % next-period value function
            for j=1:n_grid           % Solve for labor using eq. (17) in MM 
                                     % (2013)
                [n0(j,1),exitflag,n_iter,geval]=csolve('Labor_EGM',n0(j,1),[],0.000001,5,nu,gam,alpha,delta,beta,B,A,k1(j,1),z0(j,1),Wder1(j,1));
                geval_total=geval_total+geval;
            end                                                
            
            c0 = (beta*Wder1).^(-1/gam);   
                                     % Compute consumption from eq. (5) in
                                     % MM (2013)
            k0 = (B*(1-n0).^-nu./(1-alpha)/A./z0/beta./Wder1).^(1/alpha).*n0;
                                     % Compute current capital from eq. (4)
                                     % in MM (2013)
                                              
            Vder0_new = (1-delta+A*alpha*z0.*k0.^(alpha-1).*n0.^(1-alpha)).*(beta*Wder1);
                                     % Recompute the derivative of value
                                     % function in grid points

                   
            grid(:,1) = k0;          % Grid points for current capital 
             
            X0der = Polynomial_deriv_2d(grid,D); 
                                     % Construct polynomial on the current  
                                     % state variables
            vf_coef_new = X0der\Vder0_new;
                                     % Compute new coefficients for value function  
            
            feval_V_total=feval_V_total+n_grid;
 
     end % Method
        
%%%%%%%%
if Method==1 | Method==0 | Method==4
    out={V_new};
elseif Method==2 | Method==8 | Method==9
    out={vf_coef_new};
end


other_vars.k1=k1;
other_vars.c0=c0;
other_vars.n0=n0;
other_vars.k0=k0;

return
