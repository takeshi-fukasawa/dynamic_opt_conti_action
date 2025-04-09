clear
%%diary METHODSLOG.log


global n_gridk n_grida
global V k1 iter_info spectral_coef0_param lambda_param
global n0 c0
global OPI_param
global PI_linear_eq_sol_method_spec ECM_spec relative_V_spec
global analytical_EE_spec
global common_spectral_coef_spec

%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')

common_spectral_coef_spec=0;
spectral_coef0_param=1;%%% value of spectral_coef0 used in the spectral algorithm
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
spectral_coef0_param0=1;
lambda_param0=(1e-7)*(1e-3);
acceleration_spec=0;

for i = 1:6
    spectral_coef0_param=1;
    lambda_param=lambda_param0*10^(i-1);
        [out_i,other_output]=Main_function(Method,acceleration_spec,D);
        out_i=out_i(D-1,:);
        results=[results;...
            [log10(spectral_coef0_param),log10(lambda_param),out_i]];
end

%% VF-PGI-Spectral (Change lambda)
spectral_coef0_param0=1;
lambda_param0=(1e-7)*(1e-3);
acceleration_spec=1;
common_spectral_coef_spec=0;
results_spectral_change_lambda=[];

for i = 1:6
    lambda_param=lambda_param0*10^(i-1);
    spectral_coef0_param=spectral_coef0_param0;
        [out_i,other_output]=Main_function(Method,acceleration_spec,D);
        out_i=out_i(D-1,:);
        results_spectral_change_lambda=[results_spectral_change_lambda;...
            [log10(spectral_coef0_param),log10(lambda_param),out_i]];
end

%% VF-PGI-Spectral (Change spectral_coef0)
spectral_coef0_param0=1e-10;
lambda_param0=1;
acceleration_spec=1;
common_spectral_coef_spec=0;
results_spectral_change_spectral_coef0=[];

for i = 1:10
    spectral_coef0_param=spectral_coef0_param0*10^(i-1);
    lambda_param=lambda_param0;
        [out_i,other_output]=Main_function(Method,acceleration_spec,D);
        out_i=out_i(D-1,:);
        results_spectral_change_spectral_coef0=[results_spectral_change_spectral_coef0;...
            [log10(spectral_coef0_param),log10(lambda_param),out_i]];
end

%% VF-PGI-Spectral (common vs different spectral_coef)
acceleration_spec=1;
spectral_coef0_param=1;
lambda_param=1e-7;

common_spectral_coef_spec=1;
[out_common_spectral_coef,other_output]=Main_function(Method,acceleration_spec,D);
out_common_spectral_coef=out_common_spectral_coef(D-1,:);
iter_info_common_spectral_coef=iter_info;

common_spectral_coef_spec=0;
[out_different_spectral_coef,other_output]=Main_function(Method,acceleration_spec,D);
out_different_spectral_coef=out_different_spectral_coef(D-1,:);
iter_info_different_spectral_coef=iter_info;

results_spectral_coef_comparison=[out_different_spectral_coef;out_common_spectral_coef];


if save_spec==1
    filename=append('results/VF_PGI_comparison.csv');
    writematrix(round(results,3),filename)

    filename=append('results/VF_PGI_spectral_change_spectral_coef0.csv');
    writematrix(round(results_spectral_change_spectral_coef0,3),filename)

    filename=append('results/VF_PGI_spectral_change_lambda.csv');
    writematrix(round(results_spectral_change_lambda,3),filename)

    filename=append('results/VF_PGI_spectral_spectral_coef_comparison.csv');
    writematrix(round(results_spectral_coef_comparison,3),filename)
end

