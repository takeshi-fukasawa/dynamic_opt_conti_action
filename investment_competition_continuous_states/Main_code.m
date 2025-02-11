%% Steel main code
clear all

warning('off')

global elas beta_param delta_param algorithm_spec tune_param gpu_spec
global spec_precompute w_exo x_exo sd_exo
global diff lambda_param
global geval_total veval_total spec
global OPI_param
global relative_V_spec

%%% Path of spectral algorithm code
addpath('C:/Users/fukas/Dropbox/git/spectral')

addpath('./Smolyak_cpu')

%----(0) Parameters ------------------------
gpu_spec=0;
spec_precompute=1;
n_sim=1;
mu_max=3;
delta_param=0.08;%0.08;
beta_param=0.9;
elas=1.5;

resell_ratio=1;%%%%% Fixed %%%

tune_param=0;%1e-6;
n_node_exo=10;
n_node_inv=1;

n_exo=2; %% Fixed

theta=[0.12;0.3;0.12;0.3];

AR_coef=0.9*ones(1,2);
sd_exo=[0.01,0.01]*1;
sd_inv=0;%1e-4;
magnify_rate_kstk=1.5;
magnify_rate_exo=1.03;


exo_center=[4,2];
table_summary_all=[];
OPI_param=3000;

n_grid_mat=[];

for N=1:5

    d=N+2;
    k_center=ones(1,N);
    
    firm_id=N;
    
    run preliminary.m
    
    %% Simulation
    randn('seed',101)

    %% Initialize
    I_t_grid_initial0=repmat(delta_param*reshape(k_t_grid,size(k_t_grid,1),N,1)*1,1,1,n_node_inv)+...
        0.00*randn(n_grid,N,n_node_inv);
    %I_t_grid_initial0=0.05*ones(n_grid,N,n_node_inv);
    V_t_grid_initial0=V_t_grid_initial+...
        0.00*randn(n_grid,N);
    
    for relative_V_spec=0:1

        if N<=3 || relative_V_spec==1

            %% VFI algorithm
            algorithm_spec="VFI";
            run iteration.m
            
            %% VF-PGI algorithm
            algorithm_spec="gradient";
            run iteration.m
            
            %% Policy iteration
            OPI_param=3000;
            algorithm_spec="PI";
            run iteration.m
            
            if N<=3
                table_summary=round([...
                iter_results_output_func(iter_info_gradient_spectral,resid_mat_gradient_spectral);...
                iter_results_output_func(iter_info_VFI_spectral,resid_mat_VFI_spectral);...
                iter_results_output_func(iter_info_VFI,resid_mat_VFI);...
                iter_results_output_func(iter_info_PI_spectral,resid_mat_PI_spectral);...
                iter_results_output_func(iter_info_PI,resid_mat_PI);...
                 ],3);
    
            else %N>=3
                table_summary=round([...
                iter_results_output_func(iter_info_gradient_spectral,resid_mat_gradient_spectral);...
                iter_results_output_func(iter_info_VFI_spectral,resid_mat_VFI_spectral);...
                iter_results_output_func(iter_info_PI,resid_mat_PI);...
                ],3);
           end% N>=3?
    
            table_summary=[N*ones(size(table_summary,1),1),table_summary];
        
            if relative_V_spec==0
                table_summary_not_relative=table_summary;
            else
                table_summary_relative=table_summary;
            end

        else% N>=4 & relative_V_spec==0
            table_summary_not_relative=[];
        end

    end % for relative_V_spec=0:1

    table_summary_all=[table_summary_all;table_summary_not_relative;table_summary_relative];
    n_grid_mat=[n_grid_mat;[N,n_grid]];
end%N=1,2,3

writematrix(table_summary_all,append("results/results_all.csv"))
writematrix(n_grid_mat,append("results/n_grid_mat.csv"))

