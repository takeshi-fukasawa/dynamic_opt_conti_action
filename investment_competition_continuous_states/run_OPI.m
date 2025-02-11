%% OPI
table_summary_all_OPI=[];
for N=1:3
    d=N+2;
    k_center=ones(1,N);
    
    firm_id=N;
    
    run preliminary.m
    
    %% Simulation
    randn('seed',100+1)
    %% Initialize
    I_t_grid_initial0=repmat(delta_param*reshape(k_t_grid,size(k_t_grid,1),N,1)*1,1,1,n_node_inv)+...
        0.00*randn(n_grid,N,n_node_inv);
    %I_t_grid_initial0=0.05*ones(n_grid,N,n_node_inv);
    V_t_grid_initial0=V_t_grid_initial+...
        0.00*randn(n_grid,N);
    
    
    %% Policy iteration
    OPI_param=10;
   algorithm_spec="PI";
    run iteration.m   
    
    table_spectral=round([...
    iter_results_output_func(iter_info_OPI_spectral,resid_mat_OPI_spectral)],3);
    
    
    table=round([...
    iter_results_output_func(iter_info_OPI,resid_mat_OPI)],3);
    
    table_summary=round([...
    iter_results_output_func(iter_info_OPI_spectral,resid_mat_OPI_spectral);...
    iter_results_output_func(iter_info_OPI,resid_mat_OPI)],3);
    table_summary=[N*ones(size(table_summary,1),1),table_summary];
    table_summary_all_OPI=[table_summary_all_OPI;table_summary];

end%N=1,2,3

writematrix(table_summary_all_OPI,append("results/results_all_OPI.csv"))
