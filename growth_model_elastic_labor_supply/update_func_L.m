function [out,other_vars]=update_func_L(input_cell,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,spectral_spec)
  
  % Based on Maliar and Maliar (2013)
  % Modified by Takeshi Fukasawa in June 2024

   global geval_total feval_V_total n0
   global optimistic_PI_param
   global krylov_spec ECM_spec relative_V_spec analytical_EE_spec


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
   
    if 1==0 % Not precompute X1_array
       spec_V_iter=[];
       spec_V_iter.TOL=(1e-6)*max(abs(profit));
       spec_V_iter.ITER_MAX=optimistic_PI_param; 
       if spectral_spec==0
          spec_V_iter.update_spec=0;
       end
    
       if spectral_spec==3
            spec_V_iter.Anderson_acceleration=1;
       elseif spectral_spec==2
           spec_V_iter.SQUAREM_spec=1;
       end

       %%% Without precomputing X1_array
       [out,other_vars,iter_info_V_iter]=spectral_func(@VF_Bellman_L_update_func,...
           spec_V_iter,{V},...
        X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

        V_new=out{1};

    else%precompute X1_array
        X1_array=[];
        for j = 1:n_nodes                         % For each integration node...
            X1_array = cat(3,X1_array,Polynomial_2d([k1 z1(:,j)],D));   % Construct polynomial
        end
        if gam==1; uc0=log(c0); else uc0=(c0.^(1-gam)-1)/(1-gam); end
            % if gam=1, the consumption utility subfunction is
            % logarithmic; otherwise, it is power
        if nu==1; un0=log(1-n0); else un0=((1-n0).^(1-nu)-1)/(1-nu); end
            % if nu=1, the leisure utility subfunction is
            % logarithmic; otherwise, it is power
    
        profit=uc0+B*un0;
    
        V0=V;

        if krylov_spec==0
            spec_V_iter=[];
            spec_V_iter.TOL=(1e-6)*max(abs(profit));
            spec_V_iter.ITER_MAX=optimistic_PI_param; 
            
            if spectral_spec==0
                spec_V_iter.update_spec=0;
            end
            
            if spectral_spec==3
                spec_V_iter.Anderson_acceleration=1;
            elseif spectral_spec==2
                spec_V_iter.SQUAREM_spec=1;
            end

            spec_V_iter.Anderson_acceleration=1;
            [out,other_vars,iter_info_V_iter]=spectral_func(@VF_Bellman_L_given_X1_array,spec_V_iter,{V0},...
                profit,X1_array,X0,beta,weight_nodes);
            feval_V_total=feval_V_total+iter_info_V_iter.feval*n_grid;
            V_new=out{1};
        
        else%krylov_spec==1
            func_for_krylov_anonymous= @(V)func_for_krylov(V,X1_array,X0,beta,weight_nodes);
            
            ITER_MAX_gmres=100;
            TOL_gmres=1e-6;
            [V_new,flag_vec,relres,iter_gmres,resvec] = gmres(func_for_krylov_anonymous, profit,[],...
                TOL_gmres,ITER_MAX_gmres,[],[],V0); % solve for Krylov vector
            feval_V_total=feval_V_total+prod(iter_gmres)*n_grid;
        end
            
       other_vars.c0=c0;
       other_vars.n0=n0;
       other_vars.k1=k1;
       
   end
   

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
