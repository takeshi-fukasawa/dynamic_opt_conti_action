%%% Check the precision of the solution

%% Coefficients for approximation

coef_approx_V=inv_multiply_t_grid*output{2};
coef_approx_I=inv_multiply_t_grid*output{1};
[n_coef,N]=size(coef_approx_V);
n_node_inv=1;

%% Forward simulation
T=100;%%%%

k_mat=NaN(T+1,N);%t=0,...,T
exo_shock_mat=NaN(T+1,2);%t=0,...,T

resid_mat=NaN(T,2*N);%t=0,...,T

k_mat(1,:)=k_center;
exo_shock_mat(1,:)=exo_center;


random_val=randn(T,2);

for t=1:T
    k_t=k_mat(t,:);
    exo_t=exo_shock_mat(t,:);
  state_t=[k_t,exo_t];%1*(N+2)

% basis: 1*n_coef
[basis]=...
    base_func(state_t,state_min,state_max,Smol_elem,mu_max,N+2,ind);

I_t=sum(reshape(basis,n_coef,1).*reshape(coef_approx_I,n_coef,N),1);%1*N

k_mat(t+1,:)=(1-delta_param).*k_t+I_t;%1*N
exo_shock_mat(t+1,:)=AR_coef.*(exo_shock_mat(t,:)-exo_center(1,:))+exo_center(1,:)+sd_exo.*random_val(t,:);


%%% Accuracy test
V_t=sum(reshape(basis,n_coef,1).*reshape(coef_approx_V,n_coef,N),1);%1*N

[V_t1_diff,basis_t1]=V_diff_func(k_t,exo_t,coef_approx_V,...
    state_min,state_max,Smol_elem,mu_max,d,ind,w_inv);
  [inv_cost,inv_cost_diff]=inv_cost_func(state_t(1,1:N),I_t,stoch_inv_cost,theta);
    resid_I=-inv_cost_diff+beta_param*V_t1_diff;%1*N


%%%
n_pts=1;
temp=sum(reshape(basis_t1,n_pts,n_coef,n_node_inv,N).*...
    reshape(coef_approx_V,1,n_coef,1,N),2);%n_pts*1*n_node_inv*N
EV=reshape(sum(temp.*reshape(w_inv,1,1,n_node_inv,1),3),n_pts,N);%n_pts*N

% inv_cost:n_pts*N*n_node_inv
E_inv_cost=reshape(sum(inv_cost.*reshape(w_inv,1,1,n_node_inv),3),n_pts,N);

V_t_updated=pi_func(k_t,exo_t)-E_inv_cost+beta_param*EV;%n_pts*N

resid_V=V_t_updated-V_t;
resid_mat(t,:)=[resid_I,resid_V];


end % loop wrt t


