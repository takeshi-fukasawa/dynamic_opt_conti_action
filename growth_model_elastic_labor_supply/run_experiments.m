results=[];

%% VF-PGI,VF-PGI-Star,VFI,ECM,EGM,EE,PI 
krylov_spec=0;
ECM_spec=0;
for i = -2:1
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
end%i

%% EGM
i=2;
OPI_param=1;
krylov_spec=0;
[out_i,other_output]=Main_function(i,acceleration_spec,D);
results=[results;...
            [i*ones(size(out_i,1),1),out_i]];


%% EE  (Numerical/Analytical sol for solving eq)
if relative_V_spec==0
    i=3;
    for analytical_EE_spec=0:1
        [out_i,other_output]=Main_function(i,acceleration_spec,D);
        results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
    end
end

%% PI, PI-ECM

for ECM_spec=0:1
    OPI_param=3000;%sufficiently large values
    %% PI 
    krylov_spec=0;
    i=4;
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
    
    %% PI-Krylov
    krylov_spec=1;
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
    

    OPI_param=100;
    %% OPI (m=5)
    krylov_spec=0;
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
    
    %% OPI-Krylov (m=5)
    krylov_spec=1;
    OPI_param=5;
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];

    if acceleration_spec==0
        Method=6;i=6;
        [out_i,other_output]=Main_function(Method,acceleration_spec,D);
        results=[results;...
            [i*ones(size(out_i,1),1),out_i]];
    end
end

%% EGM-PI-Krylov
i=2;
krylov_spec=1;
OPI_param=3000;
[out_i,other_output]=Main_function(i,acceleration_spec,D);
results=[results;...
            [i*ones(size(out_i,1),1),out_i]];


results=results(find(results(:,2)==D),:)
