clear
%%diary METHODSLOG.log

global V k1 iter_info lambda_param alpha0_param
global common_alpha_spec

%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')


D=4;
save_spec=0;
Method=-2;

results=[];

%% Fixed point iteration and VF-PGI
alpha0_param0=1;
lambda_param0=(1e-7)*(1e-3);
spectral_spec=0;

for i = 1:6
    alpha0_param=1;
    lambda_param=lambda_param0*10^(i-1);
        [out_i,other_output]=Main_function(Method,spectral_spec,D);
        out_i=out_i(D-1,:);
        results=[results;...
            [log10(alpha0_param),log10(lambda_param),out_i]];
end

%% Spectral algorithm (Change lambda)
alpha0_param0=1;
lambda_param0=(1e-7)*(1e-3);
spectral_spec=1;
common_alpha_spec=0;
results_spectral_change_lambda=[];

for i = 1:6
    lambda_param=lambda_param0*10^(i-1);
    alpha0_param=alpha0_param0;
        [out_i,other_output]=Main_function(Method,spectral_spec,D);
        out_i=out_i(D-1,:);
        results_spectral_change_lambda=[results_spectral_change_lambda;...
            [log10(alpha0_param),log10(lambda_param),out_i]];
end

%% Spectral algorithm (Change alpha0)
alpha0_param0=1e-10;
lambda_param0=1;
spectral_spec=1;
common_alpha_spec=0;
results_spectral_change_alpha0=[];

for i = 1:10
    alpha0_param=alpha0_param0*10^(i-1);
    lambda_param=lambda_param0;
        [out_i,other_output]=Main_function(Method,spectral_spec,D);
        out_i=out_i(D-1,:);
        results_spectral_change_alpha0=[results_spectral_change_alpha0;...
            [log10(alpha0_param),log10(lambda_param),out_i]];
end

%% Spectral (common vs different alpha)
spectral_spec=1;
alpha0_param=1;
lambda_param=1e-7;

common_alpha_spec=1;
[out_common_alpha,other_output]=Main_function(Method,spectral_spec,D);
out_common_alpha=out_common_alpha(D-1,:);
iter_info_common_alpha=iter_info;

common_alpha_spec=0;
[out_different_alpha,other_output]=Main_function(Method,spectral_spec,D);
out_different_alpha=out_different_alpha(D-1,:);
iter_info_different_alpha=iter_info;

results_alpha_comparison=[out_different_alpha;out_common_alpha];


if save_spec==1
    filename=append('results/VF_PGI_comparison.csv');
    writematrix(round(results,3),filename)

    filename=append('results/VF_PGI_spectral_change_alpha0.csv');
    writematrix(round(results_spectral_change_alpha0,3),filename)

    filename=append('results/VF_PGI_spectral_change_lambda.csv');
    writematrix(round(results_spectral_change_lambda,3),filename)

    filename=append('results/VF_PGI_spectral_alpha_comparison.csv');
    writematrix(round(results_alpha_comparison,3),filename)
end

