clear
%%diary METHODSLOG.log


global n_gridk n_grida
global V k1 iter_info alpha0_param lambda_param
global n0 c0
global OPI_param
global PI_linear_eq_sol_method_spec ECM_spec relative_V_spec
global analytical_EE_spec
global common_spectral_coef_spec

%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')

common_spectral_coef_spec=0;
alpha0_param=1;%%% value of alpha0 used in the spectral algorithm
lambda_param=1e-7;%1e-12;%% Value of lambda
D=4;
n_gridk=10;n_grida=10;

OPI_param=3000;
D=4;
save_spec=1;
Method=-2;%%% VF-PGI
ECM_spec=0;
relative_V_spec=0;

results=[];

%% VF-PGI (Change lambda)
alpha0_param0=1;
lambda_param0=(1e-7)*(1e-3);
acceleration_spec=0;

for i = 1:6
    alpha0_param=1;
    lambda_param=lambda_param0*10^(i-1);
        [out_i,other_output]=Main_function(Method,acceleration_spec,D);
        out_i=out_i(D-1,:);
        results=[results;...
            [log10(alpha0_param),log10(lambda_param),out_i]];
end

%% VF-PGI-Spectral (Change lambda)
alpha0_param0=1;
lambda_param0=(1e-7)*(1e-3);
acceleration_spec=1;
common_spectral_coef_spec=0;
results_spectral_change_lambda=[];

for i = 1:6
    lambda_param=lambda_param0*10^(i-1);
    alpha0_param=alpha0_param0;
        [out_i,other_output]=Main_function(Method,acceleration_spec,D);
        out_i=out_i(D-1,:);
        results_spectral_change_lambda=[results_spectral_change_lambda;...
            [log10(alpha0_param),log10(lambda_param),out_i]];
end

%% VF-PGI-Spectral (Change alpha0)
alpha0_param0=1e-10;
lambda_param0=1;
acceleration_spec=1;
common_spectral_coef_spec=0;
results_spectral_change_alpha0=[];

for i = 1:10
    alpha0_param=alpha0_param0*10^(i-1);
    lambda_param=lambda_param0;
        [out_i,other_output]=Main_function(Method,acceleration_spec,D);
        out_i=out_i(D-1,:);
        results_spectral_change_alpha0=[results_spectral_change_alpha0;...
            [log10(alpha0_param),log10(lambda_param),out_i]];
end

%% VF-PGI-Spectral (common vs different alpha)
acceleration_spec=1;
alpha0_param=1;
lambda_param=1e-7;

common_spectral_coef_spec=1;
[out_common_alpha,other_output]=Main_function(Method,acceleration_spec,D);
out_common_alpha=out_common_alpha(D-1,:);
iter_info_common_alpha=iter_info;

common_spectral_coef_spec=0;
[out_different_alpha,other_output]=Main_function(Method,acceleration_spec,D);
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

