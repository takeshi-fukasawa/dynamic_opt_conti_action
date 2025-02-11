function out=func_for_krylov(V_t_vec,I_t,basis_t1,inv_cost,pi_mat,...
k_t,exo_t,exo_t1,basis_t,inv_multiply_t,basis_exo_t1,...
    x_inv,w_inv,...
    state_min,state_max,Smol_elem,mu_max,d,ind)

global beta_param delta_param specalgorithm_spec
global spec_precompute
global relative_V_spec

n_pts=size(I_t,1);
N=size(I_t,2);
V_t=reshape(V_t_vec,n_pts,N);

n_node_inv=size(w_inv,1);

coef_approx_V=inv_multiply_t*V_t;
[n_coef,N]=size(coef_approx_V);


% basis_t1:n_pts*n_coef*n_node_inv*N; depends on firm j !!!
% n_pts*n_coef*n_node_inv*N=> n_pts*1*n_node_inv*N

temp=sum(reshape(basis_t1,n_pts,n_coef,n_node_inv,N).*...
    reshape(coef_approx_V,1,n_coef,1,N),2);%n_pts*1*n_node_inv*N
EV=reshape(sum(temp.*reshape(w_inv,1,1,n_node_inv,1),3),n_pts,N);%n_pts*N

if relative_V_spec==1
   EV=EV-EV(1,:);
end

   out=V_t(:)-beta_param*EV(:);

end