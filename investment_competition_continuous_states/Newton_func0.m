function [I_t_updated,basis_t1,geval]=Newton_func0(...
    I_t_initial,k_t,exo_t1,basis_exo_t1,stoch_inv_cost,...
    coef_approx_V,state_min,state_max,Smol_elem,mu_max,d,ind,w_inv,theta,I_min,I_max)

global beta_param delta_param spec_precompute

global second_diff

%% Newton iteration
%I_t_initial:n_pts*N*n_node_inv; The values are used as competitiors' investment at time t

[n_pts,N,n_node_inv]=size(I_t_initial);
n_coef=size(coef_approx_V,1);

I_t_updated=I_t_initial;

basis_t1=NaN(n_pts,n_coef,n_node_inv);

ITER_MAX_Newton=30;


geval=0;
for j=1:N

    I_t_j_old=I_t_initial(:,j,:);%n_pts*1*n_node_inv
    I_t=I_t_initial;%n_pts*N*n_node_inv
    
     eps_val=1e-6;
    
    for ITER_Newton_j=1:ITER_MAX_Newton
        I_t(:,j,:)=I_t_j_old;
        I_t_eps_j=I_t;
        I_t_j_old_eps=I_t_j_old+eps_val;
        I_t_eps_j(:,j,:)=I_t_j_old_eps;
        k_t1=(1-delta_param)*k_t+I_t;%n_pts*N*n_node_inv
        k_t1_eps_j=(1-delta_param)*k_t+I_t_eps_j;%n_pts*N*n_node_inv
        
        % V_t1_diff:n_pts*N*n_node_inv
        %basis_t1:n_pts*n_coef*n_node_inv
        
        if spec_precompute==1
            [V_t1_diff_j,basis_t1_j]=...
            V_t1_diff_func_precompute(k_t1,basis_exo_t1,...
            coef_approx_V(:,j,:),state_min,state_max,Smol_elem,mu_max,ind,w_inv);

            [V_t1_diff_j_eps,basis_t1_j_eps]=...
            V_t1_diff_func_precompute(k_t1_eps_j,basis_exo_t1,...
            coef_approx_V(:,j,:),state_min,state_max,Smol_elem,mu_max,ind,w_inv);

        else
            [V_t1_diff_j,basis_t1_j]=...
            V_t1_diff_func(k_t1,exo_t1,...
            coef_approx_V(:,j,:),state_min,state_max,Smol_elem,mu_max,ind,w_inv);
        
            [V_t1_diff_j_eps,basis_t1_j_eps]=...
            V_t1_diff_func(k_t1_eps_j,exo_t1,...
            coef_approx_V(:,j,:),state_min,state_max,Smol_elem,mu_max,ind,w_inv);

        end
            %%%%%%%%%%%%%%%%%%%
            basis_t1_j=basis_t1_j(:,:,:,j);
            %%%%%%%%%%%%%%%%%%%

        [inv_cost_j,inv_cost_diff_j]=...
            inv_cost_func(k_t(:,j,:),I_t_j_old,stoch_inv_cost(:,j,:),theta);

        [inv_cost_j_eps,inv_cost_diff_j_eps]=...
            inv_cost_func(k_t(:,j,:),I_t_j_old_eps,stoch_inv_cost(:,j,:),theta);

        diff=beta_param*V_t1_diff_j(:,j,:)-inv_cost_diff_j;%n_pts*N*n_node_inv
        diff_eps=beta_param*V_t1_diff_j_eps(:,j,:)-inv_cost_diff_j_eps;%n_pts*N*n_node_inv

        second_diff=(diff_eps-diff)./eps_val;
        
        %%second_diff=min(second_diff,-0.1);%%%%
        
        scale=1.0;%%%%%
        Delta_I=-diff./second_diff.*scale;
        
        I_t_j_new=I_t_j_old+Delta_I;%n_pts*N
        
        if isempty(I_min)==0
           I_t_j_new=max(I_t_j_new,I_min);
        end
        if isempty(I_max)==0
           I_t_j_new=min(I_t_j_new,I_max);
        end
        
        DIST_Newton=max(abs(I_t_j_new(:)-I_t_j_old(:)));
        
        if ITER_Newton_j==1
            %%% Use basis_t1 based on I_t_initial
             basis_t1(:,:,:,j)=basis_t1_j;
        end
        
        if DIST_Newton<1e-10
            %DIST_Newton
            %ITER_Newton_j
            break;
        end
    
        I_t_j_old=I_t_j_new;
    
        if ITER_Newton_j==ITER_MAX_Newton
            %DIST_Newton
            %warning("Upper limit reached")
        end
    end%end of Newton loop

    geval=geval+ITER_Newton_j*n_pts*2;

    I_t_updated(:,j,:)=I_t_j_new; %%% Update I

end % for loop wrt j=1:N


if isempty(I_min)==0
   I_t_updated=max(I_t_updated,I_min);
end
if isempty(I_max)==0
   I_t_updated=min(I_t_updated,I_max);
end

return
