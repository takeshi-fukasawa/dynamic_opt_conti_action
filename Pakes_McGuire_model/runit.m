

clear all

%%% Spectral algorithm code path
addpath('C:/Users/fukas/Dropbox/git/spectral')

global diff diff_temp isentry entered

spec.ITER_MAX=300;
%spec.ITER_MAX=10;

spec.TOL=1e-6;
spec.x_min_cell={[],0};
spec.DEBUG=1;
spectral_spec=1;

%spec.common_alpha_spec=0;

if spectral_spec==0
    spec.update_spec=0;
end

% Model Parameters

c.MAX_FIRMS = 3;
c.START_FIRMS = 1;
c.EQL_TYPE = 'COMPETITION'; % COMPETITION|MONOPOLY|PLANNER

c.IND_TYPE = 'COST'; % QUALITY|COST|CAPACITY (untested, do not change)
%c.IND_TYPE = 'QUALITY'; % QUALITY|COST|CAPACITY (untested, do not change)
%c.IND_TYPE = 'CAPACITY'; % QUALITY|COST|CAPACITY (untested, do not change)

c.ENTRY_TYPE = 'RAN_ENTRY';
c.ENTRY_LOW = 0.15;
c.ENTRY_HIGH = 0.25;
c.ENTRY_SUNK = 0.2;
c.ENTRY_AT = 4;
c.BETA = 0.925;
c.BETA = 0.98;

c.DELTA = 0.7;
c.SCRAP_VAL = 0.1;
c.INV_MULT = 3;
c.INV_COST = 1;
c.QUAD_INV_COST=1.0;%%%%%%%
c.MC = 5;
c.MKT_SIZE = 5;
c.KMAX = 19;

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
profit(c);


% Solve dynamic equilibrium:
no_entry_exit_spec=0;%%%%

method=2;
[newvalue_grad,newx_grad,iter_info_grad,other_vars]=...
        eql_ma(method,no_entry_exit_spec,spec,c);

if 1==1 & c.QUAD_INV_COST==0
method=1;
[newvalue_PM,newx_PM,iter_info_PM,other_vars]=...
        eql_ma(method,no_entry_exit_spec,spec,c);
end


results_grad=[round(iter_info_grad.t_cpu,2),iter_info_grad.feval,...
    log10(iter_info_grad.DIST_table(iter_info_grad.feval,:))];

if c.QUAD_INV_COST==0
results_PM=[round(iter_info_PM.t_cpu,2),iter_info_PM.feval,...
    log10(iter_info_PM.DIST_table(iter_info_PM.feval,:))];
results=[results_PM;results_grad];

else
    results=results_grad;
end

if spectral_spec==0
    tag='';
else
    tag='_spectral';
end

if c.QUAD_INV_COST>0
    tag=append(tag,'_quad');
end

filename=append('results/PakesMcGuire_model_algorithm_comparison_summary',...
    tag,'_',string(c.MAX_FIRMS),'_',string(c.BETA),'.csv');
writematrix(round(results,2),filename)

