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
% "8": ECM-DVF
% "9": EGM-DVF

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

for relative_V_spec=0:2
    for acceleration_spec=0:1

        run run_experiments.m

        if acceleration_spec==0
            results_no_spectral=results;
        elseif acceleration_spec==1
           results_spectral=results;
        end

    end% acceleration_spec=0:1

    results=[results_spectral;results_no_spectral];

    if relative_V_spec==0
        filename=append('results/growth_model_algorithm_comparison_summary_all.csv');
    elseif relative_V_spec==1
        filename=append('results/growth_model_algorithm_comparison_summary_all_relative.csv');
    elseif relative_V_spec==2
        filename=append('results/growth_model_algorithm_comparison_summary_all_endogenous.csv');
    end

    writematrix(round(results,3),filename)
end%relative_V_spec=0:2

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SQUAREM, Anderson
relative_V_spec=0;
%%% SQUAREM
acceleration_spec=2;%SQUAREM

run run_experiments.m
results_SQUAREM=results;

%%% Anderson acceleration (One variable type)
acceleration_spec=3;%Anderson acceleration

run run_experiments.m
results_Anderson=results;

%%% Anderson acceleration (Heterogeneous variable type)
acceleration_spec=4;%Anderson acceleration

run run_experiments.m
results_Anderson_hetero=results;


filename=append('results/growth_model_algorithm_comparison_summary_alternative_acceleration.csv');

writematrix(round(...
 [results_SQUAREM;results_Anderson;results_Anderson_hetero],...
 3),filename)

%% Create summary table
results_summary=[];
for relative_V_spec=0:2
    if relative_V_spec==0
        filename=append('results/growth_model_algorithm_comparison_summary_all.csv');
    elseif relative_V_spec==1
        filename=append('results/growth_model_algorithm_comparison_summary_all_relative.csv');
    elseif relative_V_spec==2
        filename=append('results/growth_model_algorithm_comparison_summary_all_endogenous.csv');
    end
    results=readmatrix(filename);

%%%%% Check IDs!!! %%% 
    if relative_V_spec==0
        ids=[1;3;4;5;19;22;23];%VFPGI-Spectral,VFI-Spectral,VFI-ECM-Spectral,VFI-EGM-Spectral,PI,PI-Krylov
    elseif relative_V_spec>=1
        ids=[1;3;17;20;21];%VFPGI-Spectral,VFI-Spectral,PI,PI-Krylov
    end

    results_summary=[results_summary;results(ids,:)];

end%for relative_V_spec=0:2

filename=append('results/growth_model_algorithm_comparison_summary_all.csv');
results=readmatrix(filename);

ids=[15;33];% EE-Spectral, EE
results_summary=[results_summary;results(ids,:)];

filename=append('results/growth_model_algorithm_comparison_summary_concise.csv');

writematrix(results_summary,filename)

