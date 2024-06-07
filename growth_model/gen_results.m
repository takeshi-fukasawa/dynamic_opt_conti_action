clear
%%diary METHODSLOG.log

global V k1 iter_info alpha0_param lambda_param

%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')


spectral_spec=1;
common_alpha_spec=0;
alpha0_param=0.001;
lambda_param=1;
D=4;

results=[];
for i = 0:7
        [out_i,other_output]=Main_7_Methods(i,spectral_spec,D);
        results=[results;...
            [i*ones(size(out_i,1),1),out_i]];
end

results_summary=results(find(results(:,2)==D),:)

if spectral_spec==0
    tag='';
elseif spectral_spec==1
    tag='_spectral';
elseif spectral_spec==2
    tag='_SQUAREM';
end

%filename=append('results/growth_model_algorithm_comparison',...
%    tag,'.csv');

filename=append('results/growth_model_algorithm_comparison_summary',...
    tag,'.csv');

%writematrix(round(results_summary,3),filename)

%%diary off