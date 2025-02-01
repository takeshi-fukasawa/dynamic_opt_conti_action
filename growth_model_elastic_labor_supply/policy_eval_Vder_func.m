    %% Solve for value function
   [V_new,feval_V_total]=policy_eval_Vder_func(Vder0,c0,k1,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec)

global geval_total feval_V_total

    n_grid=size(V,1);

        if gam==1; uc0=log(c0); else uc0=(c0.^(1-gam)-1)/(1-gam); end
            % if gam=1, the consumption utility subfunction is
            % logarithmic; otherwise, it is power
        if nu==1; un0=log(1-n0); else un0=((1-n0).^(1-nu)-1)/(1-nu); end
            % if nu=1, the leisure utility subfunction is
            % logarithmic; otherwise, it is power
    
        profit=uc0+B*un0;

        X1der_array=[];
        for j = 1:n_nodes                         % For each integration node...
            X1der_array = cat(3,X1der_array,Polynomial_der_2d([k1 z1(:,j)],D));   % Construct polynomial
        end
    

        TOL_solve_V=1e-8;
        if krylov_spec==0
            spec_V_iter=[];
            spec_V_iter.TOL=TOL_solve_V*max(abs(profit));
            spec_V_iter.ITER_MAX=OPI_param; 
            
            spec_V_iter.update_spec=0;

            [out,other_vars,iter_info_V_iter]=spectral_func(@VF_Bellman_L_given_X1_array,spec_V_iter,{V0},...
                profit,X1_array,X0,beta,weight_nodes);
            feval_V_total=feval_V_total+iter_info_V_iter.feval*n_grid;
            V_new=out{1};
        
        else%krylov_spec==1
            func_for_krylov_anonymous= @(V)func_for_krylov(V,X1_array,X0,beta,weight_nodes);
            
            ITER_MAX_gmres=OPI_param;
            TOL_gmres=TOL_solve_V;
            [V_new,flag_vec,relres,iter_gmres,resvec] = gmres(func_for_krylov_anonymous, profit,[],...
                TOL_gmres,ITER_MAX_gmres,[],[],V0); % solve for Krylov vector
            feval_V_total=feval_V_total+prod(iter_gmres)*n_grid;
        end
       
%%%%%%%%%%
            for j = 1:n_nodes           
                X1der = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1(:,j) = X1der*vf_coef; 
                                     % Compute the derivative of value
                                     % function in the integration nodes
            end
             
            Wder1 = Vder1_EGM*weight_nodes;   
                                     % Compute expected derivative  of
                                     % next-period value function
         
            Vder0_new = (1-delta+A*alpha*z0.*k0.^(alpha-1).*n0.^(1-alpha)).*(beta*Wder1);
                                     % Recompute the derivative of value 
                                     % function in the grid points
%%%%%%%%%%%%%%%%%%

 end
   
