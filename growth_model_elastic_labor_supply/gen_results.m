clear all

global V k1 iter_info alpha0_param lambda_param
global n0 c0

%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')

spectral_spec=1;%% If 1, use the spectral algorithm. If 0, use the standard fixed point iteration.
common_alpha_spec=0;
alpha0_param=1;%%% value of alpha0 used in the spectral algorithm
lambda_param=1e-7;%% Value of lambda
D=4;

%%%%%%%%%%
results=[];
for i = -3:-2
        [out_i,other_output]=Main_function(i,spectral_spec,D);
        results=[results;...
            [i*ones(size(out_i,1),1),out_i]];
end

results_summary=results(find(results(:,2)==D),:)

if spectral_spec==0
    tag='';
elseif spectral_spec==1
    tag='_spectral';
end

filename=append('results/growth_model_algorithm_comparison_summary',...
    tag,'.csv');

writematrix(round(results_summary,3),filename)

