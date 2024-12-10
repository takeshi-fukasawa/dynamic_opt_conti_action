

clear all

addpath('C:/Users/fukas/Dropbox/git/spectral')%%%%% Spectral algorithm code path  %%
addpath('./profit_functions')
addpath('./middle_output')

global diff diff_temp isentry entered lambda_param
global wmax
global profit

spec.ITER_MAX=1000;

spec.x_min_cell={[],0};
spec.DEBUG=1;


results_all=[];
c.QUAD_INV_COST=1.0;%%%%%%% Quadratic investment cost %%%%

for N=2:2
    c.MAX_FIRMS = N; %% Number of firms %%%
    
    c.START_FIRMS = 1;
    c.EQL_TYPE = 'COMPETITION'; % COMPETITION|MONOPOLY|PLANNER
    
    %%%c.IND_TYPE = 'COST'; 
    c.IND_TYPE = 'QUALITY';%%% Consider differentiated products model as in Pakes and McGuire (1994).
    
    c.ENTRY_TYPE = 'RAN_ENTRY';
    c.ENTRY_LOW = 0.15;
    c.ENTRY_HIGH = 0.25;
    c.ENTRY_SUNK = 0.2;
    c.ENTRY_AT = 4;
    c.BETA = 0.925; %% Discount factor %%
    %%c.BETA = 0.99;
    
    c.DELTA = 0.7;
    c.SCRAP_VAL = 0.1;
    c.INV_MULT = 3;
    c.INV_COST = 1;
    c.MC = 5;
    c.MKT_SIZE = 5;
    c.KMAX = 19;
    %%c.KMAX = 20;
    
    c.WSTAR = 12;
    c.INTERCEPT = 3;
    c.FIXED_COST = 0.2;
    c.GAMMA = 1;
    c.TAU = 0.1;
    c.PROFIT_DONE = 0;
    c.EQL_DONE = 0;
    c.ACTIVE_CFG = 'default';
    
    c.PREFIX = acronym(c.EQL_TYPE, c.IND_TYPE);
    c.DS_WSTART = [c.ENTRY_AT+2; zeros(c.MAX_FIRMS-1,1)];
    c.DS_NSIMX = 10000; % 10000;
    c.DS_NSIMW = 100; % 100;
    c.DS_NRUNW = 100; % 100;
    
    
    % Compute static profit:
    %
    profit_func(c);
    
    
    % Solve dynamic equilibrium:
    no_entry_exit_spec=1;%%%% If 1, solve the model without entry/exit. If 0, solve the model with entry/exit.
    
    TOL=1e-6;
    
    %% VF-PGI-Spectral
    method="gradient";
    spec.alpha_0=1;
    lambda_param=0.01;
    TOL_vec=TOL*ones(1,2);
    TOL_vec(2)=TOL_vec(2)*lambda_param;
    spec.TOL=TOL_vec;
    
    [newvalue_grad_spectral,newx_grad_spectral,iter_info_grad_spectral,other_vars]=...
            eql_ma(method,no_entry_exit_spec,spec,c);
    
    DIST_grad_spectral=iter_info_grad_spectral.DIST_table(iter_info_grad_spectral.feval,:);
    DIST_grad_spectral(2)=DIST_grad_spectral(2)/lambda_param;
    
    if 1==0
        %% VF-PGI
        method="gradient";
        spec.alpha_0=1;
        lambda_param=0.01;
        TOL_vec=TOL*ones(1,2);
        TOL_vec(2)=TOL_vec(2)*lambda_param;
        spec.TOL=TOL_vec;
        
        spec.update_spec=0;%% Use standard fixed point iteration
        [newvalue_grad,newx_grad,iter_info_grad,other_vars]=...
                eql_ma(method,no_entry_exit_spec,spec,c);
        
        DIST_grad=iter_info_grad.DIST_table(iter_info_grad.feval,:);
        DIST_grad(2)=DIST_grad(2)/lambda_param;
    end
    
    %% PM-Spectral (VFI-star-Spectral)
    method="PM";
    spec.TOL=TOL;
    spec.alpha_0=1;
    [newvalue_PM_spectral,newx_PM_spectral,iter_info_PM_spectral,other_vars]=...
            eql_ma(method,no_entry_exit_spec,spec,c);
   
    results_grad_spectral=create_results_func(iter_info_grad_spectral);
    results_PM_spectral=create_results_func(iter_info_PM_spectral);
        
    if c.QUAD_INV_COST==0
        %% PM (VFI)
        method="PM";
        spec.TOL=TOL;
        spec.alpha_0=1;
        spec.update_spec=0; %% Use standard fixed point iteration
        [newvalue_PM,newx_PM,iter_info_PM,other_vars]=...
                eql_ma(method,no_entry_exit_spec,spec,c);
        
        results_PM=create_results_func(iter_info_PM);
    
        results=[results_grad_spectral;results_PM_spectral;results_PM];
    
    else
        results=[results_grad_spectral;results_PM_spectral];
    end
    
    
    tag=c.IND_TYPE;
    if c.QUAD_INV_COST>0
        tag=append(tag,'_quad');
    end
    results=[N*ones(size(results,1),1),results];

    filename=append('results/results_',...
        tag,'_',string(c.MAX_FIRMS),'_',string(c.BETA),'.csv');
    %writematrix(round(results,2),filename)

    results_all=[results_all;results];
end%N=1,2,3

filename=append('results/results_',...
    tag,'_',string(c.BETA),'.csv');
%writematrix(round(results_all,2),filename)

