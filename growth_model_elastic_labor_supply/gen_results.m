clear all

global V k1 iter_info alpha0_param lambda_param
global n0 c0
global OPI_param
global krylov_spec ECM_spec relative_V_spec
global analytical_EE_spec


%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')

%% Methods
% "-3": VF-PGI (update n0 and c0; treat [n0;c0] as one variable)
% "-2": VF-PGI (Updating n0 and c0)
% "-1": VF-PGI (Updating n0)
% "0" - Value function iteration
% "1" - envelope condition method iterating on value function (ECM-VF)
% "2" - endogenous grid method iterating on value function (EGM-VF)
% "3": Euler equation method (EE)
 % "4": Policy iteration method (PI) updating V
 % "5": Policy iteration method (PI) updating n0
% "6": Accelerated Value Iteration (AVI)
% "7": Safe-Accelerated Value Iteration (S-AVI)

%% Specification
% krylov_spec==0: Not apply krylov method to solve linear eq
% krylov_spec==1: Apply krylov method to solve linear eq

% relative_V_spec==0: Solve for value function
% relative_V_spec==1: Solve for relative value function
% relative_V_spec==2: Apply endogenous VFI, PI etc.

% acceleration_spec==0: Standard fixed point iteration
% acceleration_spec==1: Apply spectral algorithm
% acceleration_spec==2: Apply SQUAREM algorithm
% acceleration_spec==3: Apply Anderson acceleration (Only one variable type)
% acceleration_spec==4: Apply Anderson acceleration (Heterogeneous variable type)

% analytical_EE_spec==0: Not use analytical sol in Euler eq method
% analytical_EE_spec==1: Use analytical sol in Euler eq method

% ECM_spec==1: Use ECM mapping as the alternative to standard value
% function mapping


common_spectral_coef_spec=0;
alpha0_param=1;%%% value of alpha0 used in the spectral algorithm
lambda_param=1e-7;%1e-12;%% Value of lambda
D=4;

%% Main results
relative_V_spec=2;

for acceleration_spec=0:0

    run run_experiments.m

    if acceleration_spec==0
        results_no_spectral=results;
    elseif acceleration_spec==1
       results_spectral=results;
    elseif acceleration_spec==2
       results_Anderson=results;
    end

end% acceleration_spec

results=[results_spectral;results_no_spectral];

if relative_V_spec==0
    filename=append('results/growth_model_algorithm_comparison_summary_all.csv');
elseif relative_V_spec==1
    filename=append('results/growth_model_algorithm_comparison_summary_all_relative.csv');
elseif relative_V_spec==2
    filename=append('results/growth_model_algorithm_comparison_summary_all_endogenous.csv');
end

writematrix(round(results,3),filename)

