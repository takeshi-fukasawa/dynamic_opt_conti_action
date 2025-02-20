relative_V_spec=0;
acceleration_spec=0;
ECM_spec=0;
OPI_param=3000;
i=4;
results_PI=[];

%%% Small grid size
n_gridk=10;n_grida=10;
for PI_linear_eq_sol_method_spec=1:4
    [out,other_output]=Main_function(i,acceleration_spec,D);
    results_PI=[results_PI;...
                [n_gridk*n_grida*ones(size(out,1),1),out]];
end

%%% Large grid size
n_gridk=100;n_grida=100;
for PI_linear_eq_sol_method_spec=1:4
    [out,other_output]=Main_function(i,acceleration_spec,D);
    results_PI=[results_PI;...
                [n_gridk*n_grida*ones(size(out,1),1),out]];
end

results_PI=results_PI(find(results_PI(:,2)==D),:);
results_PI=results_PI(:,[1,4:end]);
results_PI_out=[results_PI(:,1:end-3),...
    results_PI(:,end-1)./results_PI(:,1),...% Number of iterations in policy eval step
    results_PI(:,end-1)./results_PI(:,1),...
    results_PI(:,end)];% Number of dQ_da eval

% PI-MA requires 3 matrix-vector multiplication per iteration
results_PI_out([4,8],end-1)=results_PI_out([4,8],end-1)*3;

filename=append('results/growth_model_PI_comparison.csv');

writematrix(round(results_PI_out,3),filename)



