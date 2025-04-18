clear
%%diary METHODSLOG.log

global V k1 iter_info lambda_param spectral_coef0_param
global common_spectral_coef_spec

%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')


D=4;
save_spec=0;

results=[];
results_spectral_change_spectral_coef0=[];
results_spectral_change_lambda=[];

%% Fixed point iteration and VF-PGI
spectral_coef0_param0=1;
lambda_param0=0.0000001;
spectral_spec=0;

for i = 1:6
    spectral_coef0_param=1;
    lambda_param=lambda_param0*10^(i-1);
        [out_i,other_output]=Main_7_Methods(0,spectral_spec,D);
        out_i=out_i(D-1,:);
        results=[results;...
            [log10(spectral_coef0_param),log10(lambda_param),out_i]];
end

%% Spectral algorithm (Change lambda)
spectral_coef0_param0=1;
lambda_param0=0.0000001;
spectral_spec=1;
common_spectral_coef_spec=0;

for i = 1:6
    lambda_param=lambda_param0*10^(i-1);
    spectral_coef0_param=spectral_coef0_param0;
        [out_i,other_output]=Main_7_Methods(0,spectral_spec,D);
        out_i=out_i(D-1,:);
        results_spectral_change_lambda=[results_spectral_change_lambda;...
            [log10(spectral_coef0_param),log10(lambda_param),out_i]];
end

%% Spectral algorithm (Change spectral_coef0)
spectral_coef0_param0=0.000000001;
lambda_param0=1;
spectral_spec=1;
common_spectral_coef_spec=0;

for i = 1:10
    spectral_coef0_param=spectral_coef0_param0*10^(i-1);
    lambda_param=lambda_param0;
        [out_i,other_output]=Main_7_Methods(0,spectral_spec,D);
        out_i=out_i(D-1,:);
        results_spectral_change_spectral_coef0=[results_spectral_change_spectral_coef0;...
            [log10(spectral_coef0_param),log10(lambda_param),out_i]];
end

%% Spectral (common vs different alpha)
spectral_spec=1;
spectral_coef0_param0=0.001;
lambda_param0=1;

common_spectral_coef_spec=1;
[out_common_alpha,other_output]=Main_7_Methods(0,spectral_spec,D);
out_common_alpha=out_common_alpha(D-1,:);
iter_info_common_alpha=iter_info;

common_spectral_coef_spec=0;
[out_different_alpha,other_output]=Main_7_Methods(0,spectral_spec,D);
out_different_alpha=out_different_alpha(D-1,:);
iter_info_different_alpha=iter_info;

results_alpha_comparison=[out_different_alpha;out_common_alpha];


if save_spec==1
    filename=append('results/VF_PGI_comparison.csv');
    writematrix(round(results,3),filename)

    filename=append('results/VF_PGI_spectral_change_spectral_coef0.csv');
    writematrix(round(results_spectral_change_spectral_coef0,3),filename)

    filename=append('results/VF_PGI_spectral_change_lambda.csv');
    writematrix(round(results_spectral_change_lambda,3),filename)

    filename=append('results/VF_PGI_spectral_alpha_comparison.csv');
    writematrix(round(results_alpha_comparison,3),filename)
end

