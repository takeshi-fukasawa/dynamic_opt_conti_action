%%%%%%%%%
%% Optimistic policy iteration (OPI)
OPI_param_min=3;
OPI_param_max=15;

acceleration_spec=1;
results_spectral=[];
for i=4:4
for OPI_param=OPI_param_min:OPI_param_max
           [out_i,other_output]=Main_function(i,acceleration_spec,D);
        results_spectral=[results_spectral;...
            [OPI_param*ones(size(out_i,1),1),i*ones(size(out_i,1),1),out_i]];
end%OPI_param
end%i=4:5

acceleration_spec=0;
results_no_spectral=[];
for i=4:4
for OPI_param=OPI_param_min:OPI_param_max
           [out_i,other_output]=Main_function(i,acceleration_spec,D);
        results_no_spectral=[results_no_spectral;...
            [OPI_param*ones(size(out_i,1),1),i*ones(size(out_i,1),1),out_i]];
end%OPI_param
end%i=4:5

results_OPI=[results_spectral;results_no_spectral];

results_OPI_summary=results_OPI(find(results_OPI(:,3)==D),:)

filename=append('results/growth_model_algorithm_comparison_summary_OPI_all.csv');

writematrix(round(results_OPI_summary,3),filename)

