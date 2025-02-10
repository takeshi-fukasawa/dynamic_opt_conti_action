    %% Solve for value function
function [V_new,feval_V_total]=policy_eval_func(V,n0,c0,k1,...
   X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec)

    global feval_V_total
    global PI_linear_eq_sol_method_spec OPI_param relative_V_spec

    n_grid=size(V,1);

    if gam==1; uc0=log(c0); else uc0=(c0.^(1-gam)-1)/(1-gam); end
            % if gam=1, the consumption utility subfunction is
            % logarithmic; otherwise, it is power
    if nu==1; un0=log(1-n0); else un0=((1-n0).^(1-nu)-1)/(1-nu); end
            % if nu=1, the leisure utility subfunction is
            % logarithmic; otherwise, it is power
    
    profit=uc0+B*un0;

    TOL_temp=1e-9;
    if 1==0 % Not precompute X1
       spec_V_iter=[];
       spec_V_iter.TOL=TOL_temp*max(abs(profit));
       spec_V_iter.ITER_MAX=OPI_param; 
       if acceleration_spec==0
          spec_V_iter.update_spec=0;
       end
    
       if acceleration_spec==3
            spec_V_iter.Anderson_acceleration=1;
       elseif acceleration_spec==2
           spec_V_iter.SQUAREM_spec=1;
       end

       %%% Without precomputing X1
       [out,other_vars,iter_info_V_iter]=spectral_func(@VF_Bellman_L_update_func,...
           spec_V_iter,{V},...
        X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,[],D,kdamp);

        V_new=out{1};

    else%precompute X1
        X1=[];
        for j = 1:n_nodes                         % For each integration node...
            %%X1 = cat(3,X1,Polynomial_2d([k1 z1(:,j)],D));   % Construct polynomial
            X1 = [X1;Polynomial_2d([k1 z1(:,j)],D)];   % Construct polynomial
            
        end
        %%% X1: (n_grid*n_nodes)*n_coef
    
        V0=V;

        if PI_linear_eq_sol_method_spec==0 %%% Fixed-point iteration of V operator
            spec_V_iter=[];
            spec_V_iter.TOL=TOL_temp*max(abs(profit));
            spec_V_iter.ITER_MAX=OPI_param; 
            
            spec_V_iter.update_spec=0;
            [out,other_vars,iter_info_V_iter]=spectral_func(@VF_Bellman_L_given_X1,spec_V_iter,{V0},...
                profit,X1,X0,beta,weight_nodes);
            feval_V_total=feval_V_total+iter_info_V_iter.feval*n_grid;
            V_new=out{1};
        
        elseif PI_linear_eq_sol_method_spec==1
            func_for_krylov_anonymous= @(V)func_for_krylov(V,X1,X0,beta,weight_nodes);
            
            ITER_MAX_gmres=min(OPI_param,n_grid);
            TOL_gmres=TOL_temp;

            if relative_V_spec>=1
                profit=relative_V_func(profit,relative_V_spec);
            end
            
            [V_new,flag_vec,relres,iter_gmres,resvec] = ...
                gmres(func_for_krylov_anonymous, profit,[],...
                TOL_gmres,ITER_MAX_gmres,[],[],V0); % solve for Krylov vector
            feval_V_total=feval_V_total+prod(iter_gmres)*n_grid;

         elseif PI_linear_eq_sol_method_spec>=2
   
               I=eye(n_grid);
               n_coef=size(X0,2);
               P_tilde=kron(weight_nodes',I);%n_grid*(n_grid*n_nodes)
               P_tilde_X1=sum(reshape(weight_nodes,1,n_nodes,1).*reshape(X1,n_grid,n_nodes,[]),2);%n_grid*1*n_coef
               A=I-beta*reshape(P_tilde_X1,n_grid,n_coef)*((X0'*X0)\X0');

               if PI_linear_eq_sol_method_spec==2

                    ITER_MAX_gmres=min(OPI_param,n_grid);
                    TOL_gmres=TOL_temp;
                    [V_new,flag_vec,relres,iter_gmres,resvec] = ...
                        gmres(A, profit,[],...
                        TOL_gmres,ITER_MAX_gmres,[],[],V0); % solve for Krylov vector
                    feval_V_total=feval_V_total+prod(iter_gmres)*n_grid;

               elseif PI_linear_eq_sol_method_spec==3% Exact solution of the linear equation
                    V_new=A\profit;

               elseif PI_linear_eq_sol_method_spec==4 % CG-based method proposed by Chen (2025)   
                    % A corresponds to the mapping T
                    A_trans=A';%I-t(T)

                    ITER_MAX_CG=min(OPI_param,n_grid);
                    y_i=zeros(n_grid,1);
                    s_i=profit;
                    r_i=profit;
               
                    for i=1:ITER_MAX_CG
                       alpha_i=sum(r_i.^2)/sum((A_trans*s_i).^2);
                       %%%alpha_i=sum(r_i.^2)/((A*(A_trans*s_i))'*s_i);
                       y_i_plus_1=y_i+alpha_i*s_i;
                       r_i_plus_1=profit-A*(A_trans*y_i_plus_1);
                       DIST=max(abs(r_i_plus_1));
                       if DIST<TOL_temp*max(abs(profit))
                            break;
                       end
                       beta_i=sum(r_i_plus_1.^2)/sum(r_i.^2);
                       s_i_plus_1=r_i_plus_1+beta_i*s_i;
                       
                       %%% Update variables
                       r_i=r_i_plus_1;
                       s_i=s_i_plus_1;
                       y_i=y_i_plus_1;
                    end% i=1:ITER_MAX_CG
                    V_new=A_trans*y_i;
                    
               end%PI_linear_eq_sol_method_spec==4
        end%%PI_linear_eq_sol_method_spec>=2
    end%%precompute X1 or not
return