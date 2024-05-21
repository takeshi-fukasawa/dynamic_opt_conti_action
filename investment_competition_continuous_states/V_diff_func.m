function [V_diff,basis,V_diff2]=V_diff_func(k,exo,basis_exo,coef_approx_V,...
    state_min,state_max,Smol_elem,mu_max,d,ind,w_inv)
global w_exo x_exo sd_exo
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input:
% k: n_pts*N*n_node_inv
% exo: n_pts*2
% basis_exo: n_pts*n_coef; basis at time t+1
% coef_approx_V: n_coef*N
% w_inv: n_node_inv*1

%% Output:
% V_diff: n_pts*N*n_node_inv
% basis: n_pts*n_coef*n_node_inv; basis_k.*basis_exo
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_pts,N,n_node_inv]=size(k);
n_coef=size(coef_approx_V,1);

%% Construct basis function

ns_exo=size(w_exo,2);
exo_stoch=repmat(exo,ns_exo,1)+...
    kron(sqrt(2)*sd_exo.*x_exo',ones(size(exo,1),1));%(n_pts*ns_exo)*2
state_stoch=[repmat(k,ns_exo,1),exo_stoch];%(n_pts*ns_exo)*(N+2)

% basis: n_pts*n_coef*n_node_inv
% basis_single: n_pts*n_coef*n_node_inv*(N+2)
% basis_diff: n_pts*n_coef*n_node_inv*(N+2)
% basis_diff_single: n_pts*n_coef*n_node_inv*(N+2)
if nargout<=2
[basis,basis_single,basis_diff,basis_diff_single]=...
    base_func(state_stoch,state_min,state_max,Smol_elem,mu_max,N+2,ind);
else
[basis,basis_single,basis_diff,basis_diff_single,...
    basis_diff2,basis_diff2_single]=...
    base_func(state_stoch,state_min,state_max,Smol_elem,mu_max,N+2,ind);
basis_diff2=...
    reshape(mean(reshape(basis_diff2,n_pts,ns_exo,n_coef,n_node_inv,N+2),2),...
    n_pts,n_coef,n_node_inv,N+2);%n_pts*n_coef*n_node_inv
basis_diff2_single=...
    reshape(mean(reshape(basis_diff2_single,n_pts,ns_exo,n_coef,n_node_inv,N+2),2),...
    n_pts,n_coef,n_node_inv,N+2);%n_pts*n_coef*n_node_inv*(N+2)

	end
basis=reshape(mean(reshape(basis,n_pts,ns_exo,n_coef,n_node_inv),2),...
    n_pts,n_coef,n_node_inv);%n_pts*n_coef*n_node_inv
basis_single=...
    reshape(mean(reshape(basis_single,n_pts,ns_exo,n_coef,n_node_inv,N+2),2),...
    n_pts,n_coef,n_node_inv,N+2);%n_pts*n_coef*n_node_inv*(N+2)
basis_diff=...
    reshape(mean(reshape(basis_diff,n_pts,ns_exo,n_coef,n_node_inv,N+2),2),...
    n_pts,n_coef,n_node_inv,N+2);%n_pts*n_coef*n_node_inv
basis_diff_single=...
    reshape(mean(reshape(basis_diff_single,n_pts,ns_exo,n_coef,n_node_inv,N+2),2),...
    n_pts,n_coef,n_node_inv,N+2);%n_pts*n_coef*n_node_inv*(N+2)

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
V_diff=sum(basis_diff.*reshape(coef_approx_V,1,n_coef,1,N),2);%%n_pts*1*n_node_inv*N
V_diff=permute(V_diff,[1,4,3,2]);%%n_pts*N*n_node_inv

if nargout==3
basis_diff2_adjusted=basis_diff2_single./basis_single_mean.*basis_single_mean_prod;%n_pts*n_coef*n_node_inv*N

basis_diff2=basis_diff2_adjusted(:,:,:,1:N);%n_pts*n_coef*n_node_inv*N

%%sum: n_pts*n_coef*n_node_inv*N=>n_pts*1*n_node_inv*N=>n_pts*N*n_node_inv*N
V_diff2=sum(basis_diff2.*reshape(coef_approx_V,1,n_coef,1,N),2);%%n_pts*1*n_node_inv*N
V_diff2=permute(V_diff2,[1,4,3,2]);%%n_pts*N*n_node_inv
end%if nargout==3

return
