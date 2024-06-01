function [V_diff,basis,V_diff2]=V_t1_diff_func(k_t1,exo_t1_mean,coef_approx_V,...
    state_min,state_max,Smol_elem,mu_max,ind,w_inv)
global w_exo x_exo sd_exo
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute expected E[V_t1_diff]. Expectations wrt exogenous variables

%% Input:
% k_t1: n_pts*N*n_node_inv
% exo_t1_mean: n_pts*n_exo; mean of exo_t1
% coef_approx_V: n_coef*n_var
% w_inv: n_node_inv*1

%% Output:
% V_diff: n_pts*N*n_node_inv
% basis: n_pts*n_coef*n_node_inv
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_pts,N,n_node_inv]=size(k_t1);
n_coef=size(coef_approx_V,1);
n_exo=size(exo_t1_mean,2);
n_state=N+n_exo;
n_var=size(coef_approx_V,2);

%% Construct basis function

ns_exo=size(w_exo,2);
exo_stoch=repmat(exo_t1_mean,ns_exo,1)+...
    kron(sqrt(2)*sd_exo.*x_exo',ones(size(exo_t1_mean,1),1));%(n_pts*ns_exo)*2
state_stoch=[repmat(k_t1,ns_exo,1),exo_stoch];%(n_pts*ns_exo)*(N+2)

% basis: n_pts*n_coef*n_node_inv
% basis_single: n_pts*n_coef*n_node_inv*n_state
% basis_diff: n_pts*n_coef*n_node_inv*n_state
% basis_diff_single: n_pts*n_coef*n_node_inv*n_state

if nargout<=2
    [basis_temp,basis_single_temp,basis_diff_temp,basis_diff_single_temp]=...
        base_func(state_stoch,state_min,state_max,Smol_elem,mu_max,n_state,ind);
else
    [basis_temp,basis_single_temp,basis_diff_temp,basis_diff_single_temp,...
        basis_diff2_temp,basis_diff2_single_temp]=...
        base_func(state_stoch,state_min,state_max,Smol_elem,mu_max,n_state,ind);

    %%% Sum up stochastic draws (exogenous variables)
    basis_diff2=...
        reshape(mean(reshape(basis_diff2_temp,n_pts,ns_exo,n_coef,n_node_inv,n_state),2),...
        n_pts,n_coef,n_node_inv,N+2);%n_pts*n_coef*n_node_inv
    basis_diff2_single=...
        reshape(mean(reshape(basis_diff2_single_temp,n_pts,ns_exo,n_coef,n_node_inv,n_state),2),...
        n_pts,n_coef,n_node_inv,N+2);%n_pts*n_coef*n_node_inv*n_state
end

%%% Sum up stochastic draws (exogenous variables)
basis=reshape(mean(reshape(basis_temp,n_pts,ns_exo,n_coef,n_node_inv),2),...
    n_pts,n_coef,n_node_inv);%n_pts*n_coef*n_node_inv
basis_single=...
    reshape(mean(reshape(basis_single_temp,n_pts,ns_exo,n_coef,n_node_inv,n_state),2),...
    n_pts,n_coef,n_node_inv,n_state);%n_pts*n_coef*n_node_inv*n_state
basis_diff=...
    reshape(mean(reshape(basis_diff_temp,n_pts,ns_exo,n_coef,n_node_inv,n_state),2),...
    n_pts,n_coef,n_node_inv,n_state);%n_pts*n_coef*n_node_inv
basis_diff_single=...
    reshape(mean(reshape(basis_diff_single_temp,n_pts,ns_exo,n_coef,n_node_inv,N+2),2),...
    n_pts,n_coef,n_node_inv,n_state);%n_pts*n_coef*n_node_inv*n_state

n_coef=size(basis,2);
basis_single_mean=sum(basis_single.*reshape(w_inv,1,1,n_node_inv,1),3);%n_pts*n_coef*1*N

%%% Insert 1 for zero values; needed for avoiding the division by zero
basis_single_mean(basis_single_mean==0)=1; 

basis_single_mean_prod=prod(basis_single_mean,4);%n_pts*n_coef

basis_adjusted=basis_single./basis_single_mean.*basis_single_mean_prod;%n_pts*n_coef*n_node_inv*N

basis_diff_adjusted=basis_diff_single./basis_single_mean.*basis_single_mean_prod;%n_pts*n_coef*n_node_inv*N

basis=basis_adjusted(:,:,:,1:N);%n_pts*n_coef*n_node_inv*N

basis_diff=basis_diff_adjusted(:,:,:,1:N);%n_pts*n_coef*n_node_inv*N

%%sum: n_pts*n_coef*n_node_inv*N=>n_pts*1*n_node_inv*N=>n_pts*N*n_node_inv*N
V_diff=sum(basis_diff.*reshape(coef_approx_V,1,n_coef,1,n_var),2);%%n_pts*1*n_node_inv*n_var
V_diff=permute(V_diff,[1,4,3,2]);%%n_pts*n_var*n_node_inv

if nargout==3
basis_diff2_adjusted=basis_diff2_single./basis_single_mean.*basis_single_mean_prod;%n_pts*n_coef*n_node_inv*N

basis_diff2=basis_diff2_adjusted(:,:,:,1:N);%n_pts*n_coef*n_node_inv*N

%%sum: n_pts*n_coef*n_node_inv*N=>n_pts*1*n_node_inv*N=>n_pts*N*n_node_inv*N
V_diff2=sum(basis_diff2.*reshape(coef_approx_V,1,n_coef,1,n_var),2);%%n_pts*1*n_node_inv*n_var
V_diff2=permute(V_diff2,[1,4,3,2]);%%n_pts*n_var*n_node_inv
end%if nargout==3

return
