relative_V_spec=0;
acceleration_spec=1;


ECM_spec=0;
OPI_param=3000;
i=4;
results=[];
for PI_linear_eq_sol_method_spec=1:2
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
end


