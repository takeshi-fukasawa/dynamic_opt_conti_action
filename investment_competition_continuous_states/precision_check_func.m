%%% Check the precision of the solution
function [resid_mat,k_mat,exo_shock_mat]=...
precision_check_func(I_sol,V_sol,k_center,exo_center,ind_no_precompute,AR_coef,sd_exo,theta,...
w_inv,state_min,state_max,Smol_elem,mu_max,inv_multiply_t_grid)

global beta_param delta_param

rng(100,'twister')

%% Coefficients for approximation

coef_approx_I=inv_multiply_t_grid*I_sol;
coef_approx_V=inv_multiply_t_grid*V_sol;
[n_coef,N]=size(coef_approx_V);
n_exo=size(exo_center,2);
n_state=N+n_exo;
n_node_inv=size(w_inv(:),1);

%% Forward simulation
T=100;%%%%

k_mat=NaN(T+1,N);%t=0,...,T
I_mat=NaN(T,N);%t=1,...,T

exo_shock_mat=NaN(T+1,n_exo);%t=0,...,T

resid_mat=NaN(T,2*N);%t=0,...,T

k_mat(1,:)=k_center;%%%%
exo_shock_mat(1,:)=exo_center;%%%%%


random_val=randn(T,2);

for t=1:T
    k_t=k_mat(t,:);
    exo_t=exo_shock_mat(t,:);
    state_t=[k_t,exo_t];%1*(N+2)

    % basis: 1*n_coef
    [basis_t]=...
        base_func(state_t,state_min,state_max,Smol_elem,mu_max,n_state,ind_no_precompute);
    
    I_t=sum(reshape(basis_t,n_coef,1).*reshape(coef_approx_I,n_coef,N),1);%1*N
    I_mat(t,:)=I_t;

    k_t1=(1-delta_param).*k_t+I_t;%1*N
    k_mat(t+1,:)=k_t1;

    exo_t1_mean=AR_coef.*(exo_t-exo_center)+exo_center;%1*n_exo 
    exo_shock_mat(t+1,:)=exo_t1_mean+sd_exo.*random_val(t,:);

    
    %%% Accuracy test
    V_t=sum(reshape(basis_t,n_coef,1).*reshape(coef_approx_V,n_coef,N),1);%1*N

    [basis_t1_mean,V_t1_diff]=V_t1_func(k_t1,exo_t1_mean,coef_approx_V,...
        state_min,state_max,Smol_elem,mu_max,ind_no_precompute,w_inv);
   
    stoch_inv_cost=zeros(1,N,n_node_inv);
    [inv_cost,inv_cost_diff]=inv_cost_func(k_t,I_t,stoch_inv_cost,theta);
    

    %%%
    n_pts=1;
    temp=sum(reshape(basis_t1_mean,n_pts,n_coef,n_node_inv,N).*...
        reshape(coef_approx_V,1,n_coef,1,N),2);%n_pts*1*n_node_inv*N
    EV=reshape(sum(temp.*reshape(w_inv,1,1,n_node_inv,1),3),n_pts,N);%n_pts*N

    % inv_cost:n_pts*N*n_node_inv
    E_inv_cost=reshape(sum(inv_cost.*reshape(w_inv,1,1,n_node_inv),3),n_pts,N);
    
    V_t_updated=pi_func(k_t,exo_t)-E_inv_cost+beta_param*EV;%n_pts*N

    %%% The following is not unit free??
    %resid_I=-inv_cost_diff+beta_param*V_t1_diff;%1*N
    %resid_V=V_t_updated-V_t;

    %%% Unit free
    resid_I=beta_param*V_t1_diff./inv_cost_diff-1;
    resid_V=V_t_updated./V_t-1;
    resid_mat(t,:)=[resid_I,resid_V];

end % loop wrt t


end % end function