relative_V_spec=0;
acceleration_spec=1;

krylov_spec=0;
ECM_spec=0;
i=8;
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    %results=[results;...
    %            [i*ones(size(out_i,1),1),out_i]];
    results_ECM_DVF=out_i(end,:);
k1_ECM=k1;

i=9;
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    %results=[results;...
    %            [i*ones(size(out_i,1),1),out_i]];
k1_EGM=k1;
results_EGM_DVF=out_i(end,:);
