    %% Solve for value function
function [V_new,feval_V_total]=policy_eval_func(V,n0,c0,k1,...
   X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec)

    global feval_V_total
    global krylov_spec OPI_param

    n_grid=size(V,1);

    if gam==1; uc0=log(c0); else uc0=(c0.^(1-gam)-1)/(1-gam); end
            % if gam=1, the consumption utility subfunction is
            % logarithmic; otherwise, it is power
    if nu==1; un0=log(1-n0); else un0=((1-n0).^(1-nu)-1)/(1-nu); end
            % if nu=1, the leisure utility subfunction is
            % logarithmic; otherwise, it is power
    
    profit=uc0+B*un0;

    if 1==0 % Not precompute X1_array
       spec_V_iter=[];
       spec_V_iter.TOL=(1e-6)*max(abs(profit));
       spec_V_iter.ITER_MAX=OPI_param; 
       if acceleration_spec==0
          spec_V_iter.update_spec=0;
       end
    
       if acceleration_spec==3
            spec_V_iter.Anderson_acceleration=1;
       elseif acceleration_spec==2
           spec_V_iter.SQUAREM_spec=1;
       end

       %%% Without precomputing X1_array
       [out,other_vars,iter_info_V_iter]=spectral_func(@VF_Bellman_L_update_func,...
           spec_V_iter,{V},...
        X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,[],D,kdamp);

        V_new=out{1};

    else%precompute X1_array
        X1_array=[];
        for j = 1:n_nodes                         % For each integration node...
            X1_array = cat(3,X1_array,Polynomial_2d([k1 z1(:,j)],D));   % Construct polynomial
        end
    
        V0=V;

        TOL_solve_V=1e-6;
        if krylov_spec==0
            spec_V_iter=[];
            spec_V_iter.TOL=TOL_solve_V*max(abs(profit));
            spec_V_iter.ITER_MAX=OPI_param; 
            
            spec_V_iter.update_spec=0;
            if 1==0
                if acceleration_spec==0
                    spec_V_iter.update_spec=0;
                end
                
                if acceleration_spec==3
                    spec_V_iter.Anderson_acceleration=1;
                elseif acceleration_spec==2
                    spec_V_iter.SQUAREM_spec=1;
                end
                
                spec_V_iter.Anderson_acceleration=1;

            end

            [out,other_vars,iter_info_V_iter]=spectral_func(@VF_Bellman_L_given_X1_array,spec_V_iter,{V0},...
                profit,X1_array,X0,beta,weight_nodes);
            feval_V_total=feval_V_total+iter_info_V_iter.feval*n_grid;
            V_new=out{1};
        
        else%krylov_spec==1
            func_for_krylov_anonymous= @(V)func_for_krylov(V,X1_array,X0,beta,weight_nodes);
            
            ITER_MAX_gmres=min(OPI_param,n_grid);
            TOL_gmres=TOL_solve_V;
            [V_new,flag_vec,relres,iter_gmres,resvec] = gmres(func_for_krylov_anonymous, profit,[],...
                TOL_gmres,ITER_MAX_gmres,[],[],V0); % solve for Krylov vector
            feval_V_total=feval_V_total+prod(iter_gmres)*n_grid;
        end
    end%%precompute X1_array or not
return