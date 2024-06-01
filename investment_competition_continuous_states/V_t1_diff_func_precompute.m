function [V_diff,basis,V_diff2]=V_t1_diff_func_precompute(k,basis_exo,coef_approx_V,...
    state_min,state_max,Smol_elem,mu_max,ind,w_inv)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input:
% k: n_pts*N*n_node_inv
% exo: n_pts*2
% basis_exo: n_pts*n_coef; basis at time t+1
% coef_approx_V: n_coef*n_var
% w_inv: n_node_inv*1

%% Output:
% V_diff: n_pts*n_var*n_node_inv
% basis: n_pts*n_coef*n_node_inv; basis_k.*basis_exo
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_pts,N,n_node_inv]=size(k);
n_var=size(coef_approx_V,2);


%% Construct basis function
% basis_k: n_pts*n_coef*n_node_inv
% basis_single_k: n_pts*n_coef*n_node_inv*N
% basis_diff_k: n_pts*n_coef*n_node_inv*N
% basis_diff_single_k: n_pts*n_coef*n_node_inv*N
%[basis_k,basis_single_k,basis_diff_k,basis_diff_single_k]=base_func(k,state_min(:,1:N),state_max(:,1:N),Smol_elem(:,1:N),mu_max,N);
if nargout<=2
[basis_k,basis_single_k,basis_diff_k,basis_diff_single_k]=...
    base_func(k,state_min(:,1:N),state_max(:,1:N),Smol_elem(:,1:N),mu_max,N,ind);
else
[basis_k,basis_single_k,basis_diff_k,basis_diff_single_k,...
    basis_diff2_k,basis_diff2_single_k]=...
    base_func(k,state_min(:,1:N),state_max(:,1:N),Smol_elem(:,1:N),mu_max,N,ind);
end

n_coef=size(basis_k,2);
basis_single_k_mean=sum(basis_single_k.*reshape(w_inv,1,1,n_node_inv,1),3);%n_pts*n_coef*1*N

%%% Insert 1 for zero values; needed for avoiding the division by zero
basis_single_k_mean(basis_single_k_mean==0)=1; 

basis_single_k_mean_prod=prod(basis_single_k_mean,4);%n_pts*n_coef


%% Numerical derivative
%eps=1e-5;
%[basis_k_eps,basis_single_k_eps]=base_func(k+eps,state_min(:,1:N),state_max(:,1:N),Smol_elem(:,1:N),mu_max,N);
%basis_diff_single_k_temp=(basis_single_k_eps-basis_single_k)/eps;

%% Test
%[basis_k,basis_single_k,basis_diff_k,basis_diff_single_k]=base_func(k,state_min(:,1:N),state_max(:,1:N),Smol_elem(:,1:N),mu_max,N);

%max(basis_diff_single_k_temp./basis_diff_single_k,[],'all')
%max(basis_diff_single_k_temp./basis_diff_single_k,[],'all')-...
% min(basis_diff_single_k_temp./basis_diff_single_k,[],'all')

%basis_diff_single_k=basis_diff_single_k_temp;
basis_k_adjusted=basis_single_k./basis_single_k_mean.*basis_single_k_mean_prod;%n_pts*n_coef*n_node_inv*N


basis_diff_k_adjusted=basis_diff_single_k./basis_single_k_mean.*basis_single_k_mean_prod;%n_pts*n_coef*n_node_inv*N

basis=basis_k_adjusted.*basis_exo;%n_pts*n_coef*n_node_inv*N

basis_diff=basis_diff_k_adjusted.*basis_exo;%n_pts*n_coef*n_node_inv*N

%sum(basis,'all')

%%sum:
%%n_pts*n_coef*n_node_inv*n_var=>n_pts*1*n_node_inv*n_var=>n_pts*N*n_node_inv*n_var
V_diff=sum(basis_diff.*reshape(coef_approx_V,1,n_coef,1,n_var),2);%%n_pts*1*n_node_inv*n_var
V_diff=permute(V_diff,[1,4,3,2]);%%n_pts*n_var*n_node_inv

if nargout==3
basis_diff2_k_adjusted=basis_diff2_single_k./basis_single_k_mean.*basis_single_k_mean_prod;%n_pts*n_coef*n_node_inv*N

basis_diff2=basis_diff2_k_adjusted.*basis_exo;%n_pts*n_coef*n_node_inv*N

%%sum: n_pts*n_coef*n_node_inv*N=>n_pts*1*n_node_inv*N=>n_pts*N*n_node_inv*N
V_diff2=sum(basis_diff2.*reshape(coef_approx_V,1,n_coef,1,n_var),2);%%n_pts*1*n_node_inv*n_var
V_diff2=permute(V_diff2,[1,4,3,2]);%%n_pts*n_var*n_node_inv
end%if nargout==3

return
