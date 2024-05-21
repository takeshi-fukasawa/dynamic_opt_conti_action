function pi_diff=pi_diff_func(k,exo)
%-----------------------------------------
% Purpose: differentiate pi wrt capital
%-----------------------------------------
% Input
% k:M*N;
% exo:M*2;
%-----------------------------------------
% Output
% pi_diff: M*N
%-----------------------------------------
global elas coef_firm_fixed coef_kstk gpu_spec

    [M,N]=size(k);
    demand_resid=exo(:,2)*100;
    
    if gpu_spec==1
        pi_diff=NaN(size(k,1),N,'gpuArray');%%M*N
    else
        pi_diff=NaN(size(k,1),N);%%M*N
    end
    eps=0.0001;
    for j=1:N
        k_eps=k;
        k_eps(:,j)=k_eps(:,j)+eps;
        temp=(pi_func(k_eps,exo)-pi_func(k,exo))/eps;
        pi_diff(:,j)=temp(:,j);
    end
end