%% SQUAREM
optimistic_PI_param=1000;%sufficiently large values
spectral_spec=2;%SQUAREM

results_SQUAREM=[];
for i = -3:4
    [out_i,other_output]=Main_function(i,spectral_spec,D);
    results_SQUAREM=[results_SQUAREM;...
        [i*ones(size(out_i,1),1),out_i]];
end%i

results_SQUAREM_summary=results_SQUAREM(find(results_SQUAREM(:,2)==D),:)

%% Anderson acceleration (One variable type)
optimistic_PI_param=1000;%sufficiently large values
spectral_spec=1;%Anderson acceleration

results_Anderson=[];
for i = 4:4
    [out_i,other_output]=Main_function(i,spectral_spec,D);
    results_Anderson=[results_Anderson;...
        [i*ones(size(out_i,1),1),out_i]];
end%i

results_Anderson_summary=results_Anderson(find(results_Anderson(:,2)==D),:)

%% Anderson acceleration (Heterogeneous variable type)
optimistic_PI_param=1000;%sufficiently large values
spectral_spec=4;%Anderson acceleration

results_Anderson_hetero=[];
for i = -3:4
    [out_i,other_output]=Main_function(i,spectral_spec,D);
    results_Anderson_hetero=[results_Anderson_hetero;...
        [i*ones(size(out_i,1),1),out_i]];
end%i

results_Anderson_hetero_summary=results_Anderson_hetero(find(results_Anderson_hetero(:,2)==D),:)

filename=append('results/growth_model_algorithm_comparison_summary_alternative_acceleration.csv');

writematrix(round(...
 [results_SQUAREM_summary;results_Anderson_summary;results_Anderson_hetero_summary],...
 3),filename)
