function [out,out_diff,out_diff2]=...
    inv_cost_func(kstk,inv,stoch_inv_cost,theta)

    out_plus=(theta(1)+stoch_inv_cost).*inv+theta(2)*inv.^2./kstk;
	out_minus=(theta(3)+stoch_inv_cost).*inv+theta(4)*inv.^2./kstk;
	out=out_plus.*(inv>=0)+out_minus.*(inv<0);
	
    out_diff_plus=(theta(1)+stoch_inv_cost)+2*theta(2)*inv./kstk;
    out_diff_minus=(theta(3)+stoch_inv_cost)+2*theta(4)*inv./kstk;
	out_diff=out_diff_plus.*(inv>=0)+out_diff_minus.*(inv<0);
    
    out_diff2=2*theta(2)./kstk.*(inv>0)+...
				2*theta(4)./kstk.*(inv<0);

return