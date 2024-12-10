function [out,other_vars]=...
    policy_iter_func(I_t,V_t,...
    k_t,exo_t,pi_mat,exo_t1,basis_t,inv_multiply_t,basis_exo_t1,...
    x_inv,w_inv,...
    state_min,state_max,Smol_elem,mu_max,d,ind,parameters,I_min,I_max)

% I_t: n_pts*N*n_node_inv
% V_t: n_pts*N
% k_t: n_pts*N
% exo_t1: n_pts*N*n_node_inv
% basis_t: n_pts*n_coef
% inv_multiply_t: n_pts*n_coef
% basis_exo_t1: n_pts*n_coef*1
% x_inv: n_node_inv*1
% w_inv:n_node_inv*1

%% solve optimal I (modified version)

global beta_param delta_param spec update_spec
global spec_precompute diff
global geval_total veval_total
global OPI_param

[n_pts,N,n_node_inv]=size(I_t);
n_node_inv=size(w_inv,1);

theta=parameters(1:4);
%sd_inv=parameters(end);
%stoch_inv_cost=sqrt(2)*repmat(sd_inv,1,N,1).*reshape(x_inv,1,1,n_node_inv);%1*N*n_node_inv

stoch_inv_cost=zeros(1,N,n_node_inv);

  [inv_cost,inv_cost_diff]=inv_cost_func(k_t,I_t,stoch_inv_cost,theta);


%%% Value iteration
spec_V_iter=[];
spec_V_iter.ITER_MAX=OPI_param;
spec_V_iter.TOL=1e-9;
spec_V_iter.update_spec=spec.update_spec;

[output,other_vars,iter_info]=spectral_func(...
    @V_update_func,spec_V_iter,...
    {V_t},...
I_t,inv_cost,pi_mat,...
k_t,exo_t,exo_t1,basis_t,inv_multiply_t,basis_exo_t1,...
    x_inv,w_inv,...
    state_min,state_max,Smol_elem,mu_max,d,ind);
veval_total=veval_total+n_pts*N*iter_info.n_iter;

V_t_updated=output{1};

coef_approx_V=inv_multiply_t*V_t_updated;
[I_t_updated,basis_t1,geval]=Newton_func0(I_t,k_t,exo_t1,basis_exo_t1,stoch_inv_cost,...
    coef_approx_V,state_min,state_max,Smol_elem,mu_max,d,ind,w_inv,theta,I_min,I_max);

geval_total=geval_total+geval;

out={I_t_updated,V_t_updated};
other_vars=[];

end
