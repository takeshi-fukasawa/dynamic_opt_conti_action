function basis_exo_t1=precompute_basis_exo_func(exo_t1,n_node_exo,sd_exo,Smol_elem,state_max,state_min,N,mu_max_exo,ind_exo)
global x_exo w_exo

%% Input:
% exo_t1:n_obs*2
%% Output: 
% basis_exo_t1: n_obs*n_terms

%% Purpose: Precompute basis_exo_t1

n_obs=size(exo_t1,1);
n_terms=size(Smol_elem,1);

ns_exo=size(w_exo,2);
%% generate random draws around deterministic path
exo_t1_stoch=exo_t1+sqrt(2)*reshape(sd_exo,1,2,1).*reshape(x_exo,1,2,ns_exo);

%% Construct basis_exo_t1
basis_exo_t1_stoch=base_func(exo_t1_stoch,...
    state_min(:,N+1:end),state_max(:,N+1:end),Smol_elem(:,N+1:end),mu_max_exo,2,ind_exo); %% n_obs*n_coef*n_node_exo

%% take weighted average
basis_exo_t1=reshape(sum(basis_exo_t1_stoch.*reshape(w_exo,1,1,ns_exo),3),n_obs,n_terms);

return
