function [out,other_vars]=V_update_func_given_basis_t1(V_t,I_t,basis_t1,inv_cost,pi_mat,...
k_t,exo_t,exo_t1,basis_t,inv_multiply_t,basis_exo_t1,...
    x_inv,w_inv,...
    state_min,state_max,Smol_elem,mu_max,d,ind)

global beta_param delta_param spec update_spec
global spec_precompute

[n_pts,N]=size(V_t);
n_node_inv=size(w_inv,1);

coef_approx_V=inv_multiply_t*V_t;
[n_coef,N]=size(coef_approx_V);


% basis_t1:n_pts*n_coef*n_node_inv*N; depends on firm j !!!
% n_pts*n_coef*n_node_inv*N=> n_pts*1*n_node_inv*N

temp=sum(reshape(basis_t1,n_pts,n_coef,n_node_inv,N).*...
    reshape(coef_approx_V,1,n_coef,1,N),2);%n_pts*1*n_node_inv*N
EV=reshape(sum(temp.*reshape(w_inv,1,1,n_node_inv,1),3),n_pts,N);%n_pts*N

% inv_cost:n_pts*N*n_node_inv
E_inv_cost=reshape(sum(inv_cost.*reshape(w_inv,1,1,n_node_inv),3),n_pts,N);

V_t_updated=pi_mat-E_inv_cost+beta_param*EV;%n_pts*N

if relative_V_spec==1
    V_t_updated=V_t_updated-V_t_updated(1,:);
end
 
out={V_t_updated};
other_vars=[];

end