global relative_V_spec

I_min=0.00;
I_min=[];I_max=[];

%I_min=-0.5;I_max=1.0;

spec.TOL=1e-8;
spec.spectral_coef_0=1;

    %spec=[];
    if algorithm_spec=="VFI"
        %spec.dampening_param={1.0,1.0};
    elseif algorithm_spec=="gradient"
        lambda_param=0.01;
        TOL_vec=(1e-8)*ones(1,2);
        TOL_vec(1)=TOL_vec(1)*lambda_param;
        spec.TOL=TOL_vec;
        spec.spectral_coef_0=1;
   end


I_init_val=I_t_grid_initial0;
V_init_val=V_t_grid_initial0;


%%I_init_val=I_t_gradient.*exp(0.1*rand(size(I_t_gradient)));
%%V_init_val=V_t_gradient.*exp(0.1*rand(size(I_t_gradient)));

x_max_cell={[],[]};
x_min_cell={I_min,[]};
spec.x_max_cell=x_max_cell;
spec.x_min_cell=x_min_cell;
spec.DEBUG=1;
spec.ITER_MAX=3000;

if  algorithm_spec~="PI"
    mapping=@update_func_VF;
else
    mapping=@policy_iter_func;
end

pi_mat_grid=pi_func(k_t_grid,exo_t_grid);

%%%%%%%%%% Fixed point iteration %%%%%
if N<=3 || algorithm_spec=="PI"
    geval_total=0;
    veval_total=0;
    spec.update_spec=0;
    [output,other_vars,iter_info]=spectral_func(...
        mapping,spec,...
        {I_init_val,V_init_val},...
        k_t_grid,exo_t_grid,pi_mat_grid,exo_t1_mean_grid,basis_t_grid,inv_multiply_t_grid,basis_exo_t1_mean_grid,...
       x_inv,w_inv,...
        state_min,state_max,Smol_elem,mu_max,d,ind,parameters,I_min,I_max);

    I_sol=output{1};
    V_sol=output{2};
    
    if relative_V_spec==1
        inv_cost=other_vars.inv_cost;
        [out,other_vars]=V_update_func(V_sol,I_sol,inv_cost,pi_mat_grid,...
        k_t_grid,exo_t_grid,exo_t1_mean_grid,basis_t_grid,inv_multiply_t_grid,basis_exo_t1_mean_grid,...
        x_inv,w_inv,...
        state_min,state_max,Smol_elem,mu_max,d,ind);
    
        C=(out{1}(1,:))./(1-beta_param);
        
        V_sol=V_sol-C;
    end
    
    iter_info.geval_total=geval_total;
    if veval_total==0
        iter_info.veval_total=n_grid*N*iter_info.feval;
    else
        iter_info.veval_total=veval_total;
    end
    
    [resid_mat,k_path,exo_shock_path]=...
    precision_check_func(I_sol,V_sol,k_center,exo_center,ind_no_precompute,AR_coef,sd_exo,theta,...
    w_inv,state_min,state_max,Smol_elem,mu_max,inv_multiply_t_grid);
else
    iter_info=[];I_sol=[];V_sol=[];resid_mat=[];

end % if N<=3 || algorithm_spec~="VFI"
 
%%%%%%%%%% Spectral algorithm %%%%%
if N<=3 || algorithm_spec~="PI"
    geval_total=0;
    veval_total=0;
    spec.update_spec=[];
    [output,other_vars,iter_info_spectral]=spectral_func(...
        mapping,spec,...
        {I_init_val,V_init_val},...
        k_t_grid,exo_t_grid,pi_mat_grid,exo_t1_mean_grid,basis_t_grid,inv_multiply_t_grid,basis_exo_t1_mean_grid,...
       x_inv,w_inv,...
        state_min,state_max,Smol_elem,mu_max,d,ind,parameters,I_min,I_max);
    
    I_sol_spectral=output{1};
    V_sol_spectral=output{2};
    
    if relative_V_spec==1
        inv_cost=other_vars.inv_cost;
        [out,other_vars]=V_update_func(V_sol_spectral,I_sol_spectral,inv_cost,pi_mat_grid,...
        k_t_grid,exo_t_grid,exo_t1_mean_grid,basis_t_grid,inv_multiply_t_grid,basis_exo_t1_mean_grid,...
        x_inv,w_inv,...
        state_min,state_max,Smol_elem,mu_max,d,ind);
    
        C=(out{1}(1,:))./(1-beta_param);
        V_sol_spectral=V_sol_spectral-C;
    end
    
    iter_info_spectral.geval_total=geval_total;
    if veval_total==0
        iter_info_spectral.veval_total=n_grid*N*iter_info_spectral.feval;
    else
        iter_info_spectral.veval_total=veval_total;
    end
    
    [resid_mat_spectral,k_path_spectral,exo_shock_path_spectral]=...
    precision_check_func(I_sol_spectral,V_sol_spectral,k_center,exo_center,ind_no_precompute,AR_coef,sd_exo,theta,...
    w_inv,state_min,state_max,Smol_elem,mu_max,inv_multiply_t_grid);
else
    iter_info_spectral=[];I_sol_spectral=[];V_sol_spectral=[];resid_mat_spectral=[];
end%% N<=3 || algorithm_spec~="PI"

if  algorithm_spec=="analytical"
    iter_info_analytical=iter_info;
    I_t_analytical=I_sol;
    V_t_analytical=V_sol;
    resid_mat_analytical=resid_mat;

    iter_info_analytical_spectral=iter_info_spectral;
    I_t_analytical_spectral=I_sol_spectral;
    V_t_analytical_spectral=V_sol_spectral;
    resid_mat_analytical_spectral=resid_mat_spectral;

elseif  algorithm_spec=="VFI"
    iter_info_VFI=iter_info;
    I_t_VFI=I_sol;
    V_t_VFI=V_sol;
    resid_mat_VFI=resid_mat;

    iter_info_VFI_spectral=iter_info_spectral;
    I_t_VFI_spectral=I_sol_spectral;
    V_t_VFI_spetctral=V_sol_spectral;
    resid_mat_VFI_spectral=resid_mat_spectral;
    
elseif algorithm_spec=="gradient"
    iter_info_gradient=iter_info;
    I_t_gradient=I_sol;
    V_t_gradient=V_sol;
    resid_mat_gradient=resid_mat;

    iter_info_gradient_spectral=iter_info_spectral;
    I_t_gradient_spectral=I_sol_spectral;
    V_t_gradient_spectral=V_sol_spectral;
    resid_mat_gradient_spectral=resid_mat_spectral;

elseif algorithm_spec=="PI" & OPI_param>=100
    iter_info_PI=iter_info;
    I_t_PI=I_sol;
    V_t_PI=V_sol;
    resid_mat_PI=resid_mat;

    iter_info_PI_spectral=iter_info_spectral;
    I_t_PI_spectral=I_sol_spectral;
    V_t_PI_spectral=V_sol_spectral;
    resid_mat_PI_spectral=resid_mat_spectral;

elseif algorithm_spec=="PI" & OPI_param<100
    iter_info_OPI=iter_info;
    I_t_OPI=I_sol;
    V_t_OPI=V_sol;
    resid_mat_OPI=resid_mat;

    iter_info_OPI_spectral=iter_info_spectral;
    I_t_OPI_spectral=I_sol_spectral;
    V_t_OPI_spectral=V_sol_spectral;
    resid_mat_OPI_spectral=resid_mat_spectral;
end

