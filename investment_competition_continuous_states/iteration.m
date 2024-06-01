I_min=0.00;
I_min=[];I_max=[];

%I_min=-0.5;I_max=1.0;

if 1==1
    spec=[];
    if update_spec=="PM"
        %spec.dampening_param={1.0,1.0};
    elseif update_spec=="gradient"
       spec.opt_max_spec=[0,0];
   end


   spec.alpha_0=1;

%%% {I_t_grid_initial0,V_t_grid_initial0}

I_init_val=I_t_grid_initial0;
V_init_val=V_t_grid_initial0;

if update_spec=="gradient"
    %I_init_val=I_t_analytical_spectral;
    %V_init_val=V_t_analytical_spectral;

    %I_init_val=-0.001+I_t_analytical_spectral.*exp(0.015*rand(size(I_t_analytical_spectral)));
    %V_init_val=0.1+V_t_analytical_spectral.*exp(0.01*rand(size(I_t_analytical_spectral)));

end


%%I_init_val=I_t_gradient.*exp(0.1*rand(size(I_t_gradient)));
%%V_init_val=V_t_gradient.*exp(0.1*rand(size(I_t_gradient)));

x_max_cell={[],[]};
x_min_cell={I_min,[]};
spec.x_max_cell=x_max_cell;
spec.x_min_cell=x_min_cell;
spec.DEBUG=1;
spec.ITER_MAX=500;
spec.TOL=1e-10;

%%%%%%%%%% Fixed point iteration %%%%%
spec.update_spec=0;
[output,other_vars,iter_info]=spectral_func(...
    @Bellman_I_func,spec,...
    {I_init_val,V_init_val},...
    k_t_grid,exo_t_grid,exo_t1_grid,basis_t_grid,inv_multiply_t_grid,basis_exo_t1_grid,...
   x_inv,w_inv,...
    state_min,state_max,Smol_elem,mu_max,d,ind,parameters,I_min,I_max);
end
I_sol=output{1};
V_sol=output{2};

[resid_mat,k_path,exo_shock_path]=...
precision_check_func(I_sol,V_sol,k_center,exo_center,ind_no_precompute,AR_coef,sd_exo,theta,...
w_inv,state_min,state_max,Smol_elem,mu_max,inv_multiply_t_grid);

%%%%%%%%%% Spectral algorithm %%%%%
spec.update_spec=[];
[output,other_vars,iter_info_spectral]=spectral_func(...
    @Bellman_I_func,spec,...
    {I_init_val,V_init_val},...
    k_t_grid,exo_t_grid,exo_t1_grid,basis_t_grid,inv_multiply_t_grid,basis_exo_t1_grid,...
   x_inv,w_inv,...
    state_min,state_max,Smol_elem,mu_max,d,ind,parameters,I_min,I_max);

I_sol_spectral=output{1};
V_sol_spectral=output{2};

[resid_mat_spectral,k_path_spectral,exo_shock_path_spectral]=...
precision_check_func(I_sol,V_sol,k_center,exo_center,ind_no_precompute,AR_coef,sd_exo,theta,...
w_inv,state_min,state_max,Smol_elem,mu_max,inv_multiply_t_grid);

if update_spec=="analytical"
    iter_info_analytical=iter_info;
    I_t_analytical=I_sol;
    V_t_analytical=V_sol;
    resid_mat_analytical=resid_mat;

    iter_info_analytical_spectral=iter_info_spectral;
    I_t_analytical_spectral=I_sol_spectral;
    V_t_analytical_spectral=V_sol_spectral;
    resid_mat_analytical_spectral=resid_mat_spectral;

elseif update_spec=="PM"
    iter_info_PM=iter_info;
    I_t_PM=I_sol;
    V_t_PM=V_sol;
    resid_mat_PM=resid_mat;

    iter_info_PM_spectral=iter_info_spectral;
    I_t_PM_spectral=I_sol_spectral;
    V_t_PM_spetctral=V_sol_spectral;
    resid_mat_PM_spectral=resid_mat_spectral;
    
elseif update_spec=="gradient"
    iter_info_gradient=iter_info;
    I_t_gradient=I_sol;
    V_t_gradient=V_sol;
    resid_mat_gradient=resid_mat;

    iter_info_gradient_spectral=iter_info_spectral;
    I_t_gradient_spectral=I_sol_spectral;
    V_t_gradient_spectral=V_sol_spectral;
    resid_mat_gradient_spectral=resid_mat_spectral;

end

