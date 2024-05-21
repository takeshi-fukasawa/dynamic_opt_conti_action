function [mc,mc_diff]=mc_func(k,exo,q)
%-----------------------------------------
% Purpose: Calculate mc and mc_diff given marginal cost specification
%-----------------------------------------
% Input
% k:M*N;
% exo:M*2;
%q: M*N; q==0 implies NA
%-----------------------------------------
% Output
% mc: M*N
% mc_diff: M*N
%-----------------------------------------

[M,N]=size(k);
log_mc_const=exo(:,1);

coef_kstk=-0.1;
coef_firm=0;
mc=exp(coef_kstk*log(k)+coef_firm+log_mc_const);
mc_diff=coef_kstk*mc./k;        
end