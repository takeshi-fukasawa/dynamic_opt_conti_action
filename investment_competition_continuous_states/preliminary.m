rng('default')


mu_vec=mu_max*ones(N+2,1);


%% Stochastic investment cost
[x_inv,w_inv_temp]=gausshermi(n_node_inv);
w_inv=w_inv_temp/sqrt(pi);

stoch_inv_cost=sqrt(2)*repmat(sd_inv,1,N,1).*reshape(x_inv,1,1,n_node_inv);%1*N*n_node_inv

%----------(1) generate simulated data ----------
%% Set grid point
[pi_center,mc_center,P_center,Q_center,q_center]=pi_func(k_center,exo_center);
share_center=q_center./Q_center;

k_min=k_center/magnify_rate_kstk;%1*N
k_max=k_center*magnify_rate_kstk;%1*N
exo_max=exo_center.*magnify_rate_exo;%1*2
exo_min=exo_center./magnify_rate_exo;%1*2

state_min=[k_min,exo_min];
state_max=[k_max,exo_max];

%% Construct Smolyak grid points
mu_max=max(mu_vec);
Smol_elem_iso=Smolyak_Elem_Isotrop(d,mu_max);
Smol_elem_ani=Smolyak_Elem_Anisotrop(Smol_elem_iso,mu_vec);
Smol_elem=Smol_elem_ani;
numb_terms=size(Smol_elem,1);

n_state=size(state_min,2);

%% Construct ind for state variables
row_no_precompute=reshape(repmat(1:n_state,numb_terms,1),[],1);%(numb_terms*(d=N))*1; 1,1,1,...2,2,2,...,d,d,d,...d
col_no_precompute=reshape(Smol_elem(:,1:end),[],1); %(numb_terms*N)*1; 1~m_i_max
sz_no_precompute=[N+2 numb_terms];
ind_no_precompute=sub2ind(sz_no_precompute,row_no_precompute,col_no_precompute); % (numb_terms*N)*1;

if spec_precompute==1
    %% Construct ind for kstk
    row=reshape(repmat(1:N,numb_terms,1),[],1);%(numb_terms*(d=N))*1; 1,1,1,...2,2,2,...,d,d,d,...d
    col=reshape(Smol_elem(:,1:N),[],1); %(numb_terms*N)*1; 1~m_i_max
    sz=[N numb_terms];
    ind=sub2ind(sz,row,col); % (numb_terms*N)*1;
else
    row=row_no_precompute;
    col=col_no_precompute;
    sz=sz_no_precompute;
    ind=ind_no_precompute;
end

%% Construct ind for exo
row=reshape(repmat(1:2,numb_terms,1),[],1);%(numb_terms*(d=2))*1; 1,1,1,...2,2,2,...,d,d,d,...d
col=reshape(Smol_elem(:,N+1:end),[],1); %(numb_terms*(d=2))*1; 1~m_i_max
sz=[2 numb_terms];
ind_exo=sub2ind(sz,row,col); % (numb_terms*2)*1;

Smol_grid=Smolyak_Grid(d,mu_vec,Smol_elem);

k_grid=Rescale(Smol_grid(:,1:N),k_min,k_max,-1,1);%M*N
exo_grid=Rescale(Smol_grid(:,(N+1):end),exo_min,exo_max,-1,1);%M*2

[pi_grid,mc_grid,P_grid,Q_grid,q_grid]=pi_func(k_grid,exo_grid);
share_grid=q_grid./Q_grid;

%--------  Set theta1 -----------------------%
%%% Assume that theta is approximately consistent with k_center.
k_t1_desired=k_center;
pi_t1_diff_initial=pi_diff_func(k_t1_desired,exo_center);%%1*N

theta_temp=...
	(beta_param*pi_t1_diff_initial(1)+beta_param*theta(2)*delta_param^2)/(1-beta_param*(1-delta_param))...
	-2*theta(2)*delta_param;
theta(1)=mean(theta_temp);

theta(3)=theta(1)*resell_ratio;%0.99;
theta(4)=theta(2);

%----- Set initial values of variables for value function approximation ----------------%
basis_t_grid=Smolyak_Polynomial(Smol_grid,d,mu_max,Smol_elem,[]);%M*n_coef

[n_grid,n_coef]=size(basis_t_grid);

%%% Use terminal state condition as initial values of V (Assume that capital stocks not change from the next period)
inv_cost=inv_cost_func(k_grid,k_grid*delta_param,stoch_inv_cost,theta);
inv_cost_mean=sum(inv_cost.*reshape(w_inv,1,1,n_node_inv),3);
V_t_grid_initial=(pi_func(k_grid,exo_grid)-inv_cost_mean)/(1-beta_param);
V_t_grid_initial_terminal=V_t_grid_initial;

%% Calculate the coefficient for approximation
inv_multiply_t_grid=inv(basis_t_grid'*basis_t_grid+tune_param*eye(n_coef,n_coef))*basis_t_grid';

coef_approx_V_initial=inv_multiply_t_grid*V_t_grid_initial; %% M*N

corr_check_initial=diag(corr(V_t_grid_initial,basis_t_grid*coef_approx_V_initial))

%% Validation of the code
%run validation.m

k_t_grid=k_grid;%%M*N
exo_t_grid=exo_grid;%%M*2
exo_t1_grid=AR_coef.*(exo_t_grid-exo_center)+exo_center;%%M*2

%% Precompute basis_exo_t1
[x_temp,w_temp]=gausshermi(n_node_exo);
w=w_temp/sqrt(pi);

[X_temp,Y_temp]=meshgrid(x_temp,x_temp);
x_exo=[X_temp(:),Y_temp(:)]'*sqrt(2);%2*(n_node_exo^2==ns)

[X_temp,Y_temp]=meshgrid(w,w);
weight_mat=[X_temp(:),Y_temp(:)]';%2*(n_node^2==ns)
weight=prod(weight_mat,1);%1*ns
w_exo=weight./sum(weight,2);%1*ns;necessary??


basis_exo_t1_grid=precompute_basis_exo_func(exo_t1_grid,n_node_exo,sd_exo,Smol_elem,state_max,state_min,N,mu_max,ind_exo);

parameters=[theta;sd_inv];
