clear
%%diary METHODSLOG.log

global V k1 iter_info lambda_param alpha0_param
global common_alpha_spec

addpath('C:/Users/fukas/Dropbox/git/spectral')


D=4;

results=[];
results_spectral=[];

%% Fixed point iteration and VF-PGI
alpha0_param0=1;
lambda_param0=0.0000001;
spectral_spec=0;

for i = 1:6
    alpha0_param=1;
    lambda_param=lambda_param0*10^(i-1);
        [out_i,other_output]=Main_7_Methods(0,spectral_spec,D);
        out_i=out_i(D-1,:);
        results=[results;...
            [spectral_spec,alpha0_param,lambda_param,out_i]];
end

%% Spectral algorithm
alpha0_param0=0.000000001;
lambda_param0=1;
spectral_spec=1;
common_alpha_spec=0;

for i = 1:10
    alpha0_param=alpha0_param0*10^(i-1);
    lambda_param=lambda_param0;
        [out_i,other_output]=Main_7_Methods(0,spectral_spec,D);
        out_i=out_i(D-1,:);
        results_spectral=[results_spectral;...
            [spectral_spec,alpha0_param,lambda_param,out_i]];
end

%% Spectral (common vs different alpha)
spectral_spec=1;
alpha0_param0=0.001;
lambda_param0=1;

common_alpha_spec=0;
[out_common_alpha,other_output]=Main_7_Methods(0,spectral_spec,D);
out_common_alpha=out_common_alpha(D-1,:);

common_alpha_spec=1;
lambda_param0=0.01;
[out_different_alpha,other_output]=Main_7_Methods(0,spectral_spec,D);
out_different_alpha=out_different_alpha(D-1,:);
iter_info_common_alpha=iter_info;

results_alpha_comparison=[out_different_alpha;out_common_alpha];

if spectral_spec==0
    tag='';
elseif spectral_spec==1
    tag='_spectral';
elseif spectral_spec==2
    tag='_SQUAREM';
end

%filename=append('growth_model_algorithm_comparison',...
%    tag,'.csv');

filename=append('growth_model_algorithm_comparison_summary',...
    tag,'.csv');

%writematrix(round(results_summary,3),filename)
