function [Q_diff]=Q_diff_func(...
    I_t_j_state,state_id,firm_id,I_t,k_t,exo_t1,basis_exo_t1,stoch_inv_cost,...
    coef_approx_V,state_min,state_max,Smol_elem,mu_max,d,ind,w_inv,theta,I_min,I_max)

    %%% Not validated
    
    I_t_state=I_t(state_id,:,:);
    I_t_state(:,firm_id,:)=I_t_initial_j;

    k_t_state=k_t(state_id,:,:);
    k_t1_state=(1-delta_param)*k_t_state+I_t_state;%n_pts*N*n_node_inv

    if spec_precompute==1
        [basis_t1_j,V_t1_diff_j]=...
        V_t1_func_precompute(k_t1,basis_exo_t1,...
        coef_approx_V(:,j,:),state_min,state_max,Smol_elem,mu_max,ind,w_inv);
    else
        [basis_t1_j,V_t1_diff_j]=...
        V_t1_func(k_t1,exo_t1,...
        coef_approx_V(:,j,:),state_min,state_max,Smol_elem,mu_max,ind,w_inv);
    end

    basis_t1_j=basis_t1_j(:,:,:,j);

    [inv_cost_j,inv_cost_diff_j]=...
        inv_cost_func(k_t(:,j,:),I_t_j_old,stoch_inv_cost(:,j,:),theta);

    Q_diff=beta_param*V_t1_diff_j(:,j,:)-inv_cost_diff_j;%n_pts*N*n_node_inv


end
