clear all

global V k1 iter_info alpha0_param lambda_param
global n0 c0
global optimistic_PI_param

%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')

common_alpha_spec=0;
alpha0_param=1;%%% value of alpha0 used in the spectral algorithm
lambda_param=1e-7;%% Value of lambda
D=4;

%%%%%%%%%%
%% Main results
optimistic_PI_param=500;%sufficiently large values

for spectral_spec=1:1
results=[];
for i = 5:5
        [out_i,other_output]=Main_function(i,spectral_spec,D);
        results=[results;...
            [i*ones(size(out_i,1),1),out_i]];
end%i
if spectral_spec==0
   results_no_spectral=results;
else
   results_spectral=results;
end
end% spectral_spec
results=[results_spectral;results_no_spectral];

results_summary=results(find(results(:,2)==D),:)

filename=append('results/growth_model_algorithm_comparison_summary_all.csv');

%writematrix(round(results_summary,3),filename)

%%%%%%%%%
%% Optimistic policy iteration (OPI)
spectral_spec=1;
results_spectral=[];
for i=4:5
for optimistic_PI_param=1:5
           [out_i,other_output]=Main_function(i,spectral_spec,D);
        results_spectral=[results_spectral;...
            [i*ones(size(out_i,1),1),out_i]];
end%optimistic_PI_param
end%i=4:5

spectral_spec=0;
results_no_spectral=[];
for i=4:5
for optimistic_PI_param=1:5
           [out_i,other_output]=Main_function(i,spectral_spec,D);
        results_no_spectral=[results_no_spectral;...
            [i*ones(size(out_i,1),1),out_i]];
end%optimistic_PI_param
end%i=4:5

results_OPI=[results_spectral;results_no_spectral];

results_OPI_summary=results(find(results_OPI(:,2)==D),:)

filename=append('results/growth_model_algorithm_comparison_summary_OPI_all.csv');

%writematrix(round(results_OPI_summary,3),filename)
