results=[];

%if 1==0%%%%%%%%%%%%%%%
%% VF-PGI,VF-PGI-Star,VFI
krylov_spec=0;
ECM_spec=0;
OPI_param=1;
for i = -2:1
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
end%i

%% VFI-EGM
i=2;
krylov_spec=0;
[out_i,other_output]=Main_function(i,acceleration_spec,D);
results=[results;...
            [i*ones(size(out_i,1),1),out_i]];

%end%%%%%%%%%%%


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
    

    %% OPI
    OPI_param=100;
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
end% for ECM_spec=0:1

%% PI-Krylov-EGM
i=2;
krylov_spec=1;
OPI_param=3000;
[out_i,other_output]=Main_function(i,acceleration_spec,D);
results=[results;...
            [i*ones(size(out_i,1),1),out_i]];

%% AVFI, AVFI-ECM
if acceleration_spec==0
    for ECM_spec=0:1
        Method=6;i=6;
        [out_i,other_output]=Main_function(Method,acceleration_spec,D);
        results=[results;...
            [i*ones(size(out_i,1),1),out_i]];
    end%ECM_spec=0:1
end%if acceleration_spec==0

%% DVFI-ECM, DVFI-EGM
if relative_V_spec==0
   for i=8:9
              [out_i,other_output]=Main_function(i,acceleration_spec,D);
        results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
    end
end

%% EE  (Numerical/Analytical sol for solving eq)
if relative_V_spec==0 || acceleration_spec>=2
    i=3;
    for analytical_EE_spec=0:1
        [out_i,other_output]=Main_function(i,acceleration_spec,D);
        results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
    end
end


results=results(find(results(:,2)==D),:)
