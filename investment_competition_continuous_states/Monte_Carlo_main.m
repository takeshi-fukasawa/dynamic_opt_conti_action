%% Steel main code
clear all

global elas beta_param delta_param update_spec tune_param gpu_spec
global spec_precompute w_exo x_exo sd_exo
global diff lambda_param

DEBUG=0;
addpath('./Smolyak_cpu')

%%% Path of spectral algorithm code
addpath('C:/Users/fukas/Dropbox/git/spectral')

%----(0) Parameters ------------------------
lambda_param=1;

gpu_spec=0;
spec_precompute=1;
n_sim=1;
N=2;d=N+2;
mu_max=2;
delta_param=0.08;%0.08;
beta_param=0.9;
elas=1.5;

%%resell_ratio=0.8;%1;
resell_ratio=1;%%%%%

tune_param=0;%1e-6;
n_node_exo=5;
n_node_inv=1;

n_exo=2; %% Fixed

adj_cost=0.3;
theta=[0.2;adj_cost;0.2;adj_cost];

AR_coef=0.9*ones(1,2);
sd_exo=[0.01,0.01]*1;
sd_inv=0;%1e-4;
magnify_rate_kstk=1.8;
magnify_rate_exo=1.1;
k_center=ones(1,N);
exo_center=[4,2];%%%%


firm_id=N;

run preliminary.m

%% Simulation
for i=1:n_sim
randn('seed',100+i)
%% Initialize
I_t_grid_initial0=repmat(delta_param*reshape(k_t_grid,size(k_t_grid,1),N,1)*1,1,1,n_node_inv)+...
    0.01*randn(n_grid,N,n_node_inv);
%I_t_grid_initial0=0.05*ones(n_grid,N,n_node_inv);
V_t_grid_initial0=V_t_grid_initial+...
    0.00*randn(n_grid,N);
%V_t_grid_initial0=0.001*randn(n_grid,N);
%V_t_grid_initial=randn(n_grid,N);
%I_t_grid_initial=repmat(delta_param*reshape(k_t_grid,size(k_t_grid,1),N,1)*1,1,1,n_node_inv)+...
%    0.1*randn(n_grid,N,n_node_inv);

%% Pakes McGuire(1994) algorithm
if 1==1
update_spec="PM";
tic
run iteration.m
end%%%%%


%% New algorithm
update_spec="analytical";
tic
run iteration.m


%% Gradient algorithm
update_spec="gradient";
tic
run iteration.m

end


%max(V_t_grid_initial0./V_t_analytical_algorithm)
%max(I_t_grid_initial0-I_t_analytical_algorithm)
max(I_t_PM./I_t_analytical)
max(V_t_PM./V_t_analytical)


table_spectral=round([...
iter_results_output_func(iter_info_gradient_spectral,resid_mat_gradient_spectral);...
iter_results_output_func(iter_info_PM_spectral,resid_mat_PM_spectral);...
iter_results_output_func(iter_info_analytical_spectral,resid_mat_analytical_spectral)],3);


table=round([...
iter_results_output_func(iter_info_gradient,resid_mat_gradient);...
iter_results_output_func(iter_info_PM,resid_mat_PM);...
iter_results_output_func(iter_info_analytical,resid_mat_analytical)],3);



%writematrix([table;table_spectral],append("results/results_",string(N),".csv"))

