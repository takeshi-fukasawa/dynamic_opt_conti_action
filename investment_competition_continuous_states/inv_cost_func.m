function [out,out_diff,out_diff2]=...
    inv_cost_func(kstk,inv,stoch_inv_cost,theta)

    out=(theta(1)+stoch_inv_cost).*inv+theta(2)*inv.^2./kstk;
	
    out_diff=(theta(1)+stoch_inv_cost)+2*theta(2)*inv./kstk;
    
    out_diff2=2*theta(2)./kstk;

return