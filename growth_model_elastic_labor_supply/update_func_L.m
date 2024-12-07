function [out,other_vars]=update_func_L(input_cell,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,spectral_spec)
  
  % Based on Maliar and Maliar (2013)
  % Modified by Takeshi Fukasawa in June 2024

   global geval_total feval_V_total n0
   global optimistic_PI_param

if Method==1 | Method==0 | Method==4
    V=input_cell;
elseif Method==2
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
                            
        [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                % Recompute value function using the Bellman 
                % equation
        
         grid(:,1) = k0;        % Grid points for current capital 
            
         X0 = Polynomial_2d(grid,D);   % Construct polynomial on 
                                          % current state variables
         vf_coef_new=X0\V_new;

         vf_coef_new = kdamp*vf_coef_new + (1-kdamp)*vf_coef;   
                % Update Vvf_coef using damping

elseif Method==4
     
%==================================================================
% Method 4. Policy iteration (PI) directly updating V
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

   spec_V_iter=[];
   spec_V_iter.TOL=1e-6;
   spec.V_iter.ITER_MAX=optimistic_PI_param; 
     if spectral_spec==0
      spec.update_spec=0;
    end
 [out,other_vars,iter_info_V_iter]=spectral_func(@V_update_func,spec_V_iter,{V},...
    X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

   feval_V_total=feval_V_total+iter_info_V_iter.feval*n_grid;

   V_new=out{1};

end
        
%%%%%%%%
if Method==1 | Method==0 | Method==4
    out={V_new};
elseif Method==2
    out={vf_coef_new};
end


other_vars.k1=k1;
other_vars.c0=c0;
other_vars.n0=n0;
other_vars.k0=k0;

return
