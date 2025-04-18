function [out,other_vars]=...
    update_func_VF(I_t,V_t,...
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

global beta_param delta_param spec algorithm_spec
global spec_precompute diff
global geval_total
global relative_V_spec

[n_pts,N,n_node_inv]=size(I_t);
n_node_inv=size(w_inv,1);

coef_approx_V=inv_multiply_t*V_t;
[n_coef,N]=size(coef_approx_V);

theta=parameters(1:4);
%sd_inv=parameters(end);
%stoch_inv_cost=sqrt(2)*repmat(sd_inv,1,N,1).*reshape(x_inv,1,1,n_node_inv);%1*N*n_node_inv

stoch_inv_cost=zeros(1,N,n_node_inv);


%% --------------- (1) Calculate FOC (investment choice problem) -----------------
if algorithm_spec=="analytical"|algorithm_spec=="gradient"
    % V_t1_diff:n_pts*N*n_node_inv
    if spec_precompute==1
        [basis_t1,V_t1_diff]=...
            V_t1_func_precompute((1-delta_param)*reshape(k_t,n_pts,N,1)+I_t,...
            basis_exo_t1,coef_approx_V,state_min,state_max,Smol_elem,mu_max,ind,w_inv);

    else
         [basis_t1,V_t1_diff]=...
            V_t1_func((1-delta_param)*reshape(k_t,n_pts,N,1)+I_t,...
        exo_t1,coef_approx_V,state_min,state_max,Smol_elem,mu_max,ind,w_inv);
    end

    [inv_cost,inv_cost_diff]=inv_cost_func(k_t,I_t,stoch_inv_cost,theta);
    I_t_updated=I_update_func(I_t,k_t,V_t1_diff,stoch_inv_cost,theta,algorithm_spec,I_min,I_max);
    
    %%% The following is unstable????
    %%%[inv_cost,inv_cost_diff]=inv_cost_func(k_t,I_t_updated,stoch_inv_cost,theta);

    if algorithm_spec=="analytical"
        %diff=-inv_cost_diff+beta_param*V_t1_diff;
    end

    geval_total=geval_total+n_pts*N;

else %% VFI algorithm
    %%% basis_t1: based on I_t, rather than I_t_updated
    [I_t_updated,basis_t1,geval]=Newton_func0(I_t,k_t,exo_t1,basis_exo_t1,stoch_inv_cost,...
    coef_approx_V,state_min,state_max,Smol_elem,mu_max,d,ind,w_inv,theta,I_min,I_max);
  
    geval_total=geval_total+geval;

    %%%[inv_cost,inv_cost_diff]=inv_cost_func(k_t,I_t,stoch_inv_cost,theta);
    [inv_cost,inv_cost_diff]=inv_cost_func(k_t,I_t_updated,stoch_inv_cost,theta);
end

% basis_t1:n_pts*n_coef*n_node_inv*N; depends on firm j !!!
% n_pts*n_coef*n_node_inv*N=> n_pts*1*n_node_inv*N

temp=sum(reshape(basis_t1,n_pts,n_coef,n_node_inv,N).*...
    reshape(coef_approx_V,1,n_coef,1,N),2);%n_pts*1*n_node_inv*N
EV=reshape(sum(temp.*reshape(w_inv,1,1,n_node_inv,1),3),n_pts,N);%n_pts*N

% inv_cost:n_pts*N*n_node_inv
E_inv_cost=reshape(sum(inv_cost.*reshape(w_inv,1,1,n_node_inv),3),n_pts,N);

V_t_updated=pi_mat-E_inv_cost+beta_param*EV;%n_pts*N

if relative_V_spec==1
    V_t_updated=V_t_updated-V_t_updated(1,:);%%% Relative value function iteration (cf. Bray 2019)
end

out={I_t_updated,V_t_updated};

other_vars=[];
other_vars.inv_cost=inv_cost;

return