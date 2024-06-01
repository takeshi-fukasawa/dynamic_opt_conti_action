function [pi,mc,P,Q,q]=pi_func(k_grid,exo_grid)
%-----------------------------------------
% Purpose: Calculate profit pi
%-----------------------------------------
% Input
% k_grid:M*N;
% exo_grid:M*2;
%-----------------------------------------
% Output
% pi: M*N
% mc: M*N
% P: M*1
% Q: M*1
% q: M*N
%-----------------------------------------
global elas
[M,N]=size(k_grid);


log_demand_resid=exo_grid(:,2);

mc=mc_func(k_grid,exo_grid,0);
P=sum(mc,2)/(N-1/elas);
Q=exp(log_demand_resid).*P.^(-elas);
q=(P-mc).*elas.*Q./P;
pi=(P-mc).*q;

end
