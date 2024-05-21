%% Validation of the main code

%% Validating the accuracy of approximation
%% Cross validation 
outside_ratio_k=0.5;
outside_ratio_exo=0.5;
M_valid=100;
rng(100);k_valid=k_min+(k_max-k_min).*rand(M_valid,N);
rng(200);exo_valid=exo_min+outside_ratio_exo*(exo_max-exo_min).*rand(M_valid,2);
basis_valid=base_func([k_valid,exo_valid],state_min,state_max,Smol_elem,mu_max,d,[]);%%M*n_basis;

inv_cost_mean=sum(inv_cost_func(k_valid,k_valid*delta_param,stoch_inv_cost,theta).*reshape(w_inv,1,1,n_node_inv),3);
V_valid=(pi_func(k_valid,exo_valid)-inv_cost_mean)/(1-beta_param);
V_valid_predict=basis_valid*coef_approx_V_initial;
corr_approx_cros_valid=diag(corr(V_valid(:,1),V_valid_predict(:,1))) %% should be close to 1

%hist(V_valid./V_valid_predict)

%% check the validity of analytical derivative of pi
eps=1e-4;
k_grid_eps=k_grid;
k_grid_eps(:,1)=k_grid(:,1)+eps;
pi_t1=pi_func(k_grid,exo_grid);
pi_t1_eps=pi_func(k_grid_eps,exo_grid);
diff_pi_temp=(pi_t1_eps-pi_t1)/eps;
pi_t1_diff=pi_diff_func(k_grid,exo_grid);
check_pi_derivative=pi_t1_diff(:,1)./diff_pi_temp(:,1);%% should be equal to one, except for both zero case (NaN)


%% second derivative of pi func
eps=0.00001;
%k_t1_desired=k_t1_desired*0.5;
%k_t1_desired_eps=k_t1_desired;
%k_t1_desired_eps(firm_id)=k_t1_desired_eps(firm_id)+eps;
%pi_t1_diff_initial=pi_diff_func(k_t1_desired,exo_center);%%1*N
%pi_t1_diff_eps=pi_diff_func(k_t1_desired_eps,exo_center);%%1*N
%second_derivative_pi=(pi_t1_diff_eps(firm_id)-pi_t1_diff_initial(firm_id))/eps %% strictly speaking, not SOC if theta(2,3)!=0
%second_derivative_Euler=(second_derivative_pi-2*theta(2)*delta_param^2)/(1-beta_param)...
%    -2*theta(2)-2*theta(3)./k_t1_desired

%[mc_temp,mc_diff_temp]=mc_func(k_t1_desired,exo_center);
%[mc_temp2,mc_diff_temp2]=mc_func(k_t1_desired_eps,exo_center);
%mc_diff2=(mc_diff_temp2-mc_diff_temp)/eps;%% second derivative of mc

%% Shape of pi_func
k_x=[1:20]';
mc_y=mc_func(k_x,exo_center,0);
%plot(k_x,mc_y)
pi_y=pi_func([k_x,10*ones(20,1)],exo_center);
%plot(k_x,pi_y)


%% Validate basis_diff
[basis_k,basis_single_k,...
    basis_diff_k,basis_diff_single_k,...
    basis_diff2_k,basis_diff2_single_k]=...
    base_func(k_grid,state_min(:,1:N),state_max(:,1:N),...
    Smol_elem(:,1:N),mu_max,N,ind);

eps=1e-6;

for j=1:N
k_eps=k_grid;
k_eps(:,j)=k_eps(:,j)+eps;
[basis_k_eps,temp,basis_diff_k_eps]=base_func(k_eps,state_min(:,1:N),state_max(:,1:N),...
    Smol_elem(:,1:N),mu_max,N,ind);%n_pts*n_terms
basis_diff_k_numer(:,:,j)=(basis_k_eps-basis_k)/eps;
basis_diff2_k_numer(:,:,:,:,j)=(basis_diff_k_eps-basis_diff_k)/eps;

end

%%% If diff_analytical <1e^6, assign zero/zero
diff_analytical=gather(squeeze(basis_diff_k(:,:,:,1)));
diff_numer=gather(basis_diff_k_numer(:,:,1));
ratio=diff_numer./diff_analytical.*(diff_analytical>1e-6);
max(diff_analytical-diff_numer)./max(diff_analytical)%%Should be sufficiently close tp zero

