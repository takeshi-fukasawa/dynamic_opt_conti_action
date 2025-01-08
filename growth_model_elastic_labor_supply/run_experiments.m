results=[];

%% VF-PGI,VF-PGI-Star,VFI,ECM,EGM,EE,PI 
krylov_spec=0;
ECM_spec=0;
for i = -2:3
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
end%i

%% PI, PI-ECM

for ECM_spec=0:1
    optimistic_PI_param=3000;%sufficiently large values
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
    

    optimistic_PI_param=5;
    %% OPI (m=5)
    krylov_spec=0;
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
    
    %% OPI-Krylov (m=5)
    krylov_spec=1;
    optimistic_PI_param=5;
    [out_i,other_output]=Main_function(i,acceleration_spec,D);
    results=[results;...
                [i*ones(size(out_i,1),1),out_i]];
end

if acceleration_spec==0
   Method=6;i=6;
   [out_i,other_output]=Main_function(Method,acceleration_spec,D);
   results=[results;...
         [i*ones(size(out_i,1),1),out_i]];
end

results=results(find(results(:,2)==D),:)
