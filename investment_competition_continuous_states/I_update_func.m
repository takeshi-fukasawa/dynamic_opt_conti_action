function I_updated=I_update_func(I_initial,k_t,V_t1_diff,stoch_inv_cost,theta,update_spec,I_min,I_max)
global beta_param lambda_param

if update_spec=="analytical"
    I_opt_plus=(beta_param*V_t1_diff-theta(1)-stoch_inv_cost)./(2*theta(2)./k_t);
    I_opt_minus=(beta_param*V_t1_diff-theta(3)-stoch_inv_cost)./(2*theta(4)./k_t);
    plus_dummy=(beta_param*V_t1_diff>theta(1)+stoch_inv_cost);
    minus_dummy=(beta_param*V_t1_diff<theta(3)+stoch_inv_cost);
    zero_dummy=1-plus_dummy-minus_dummy;
    I_updated=I_opt_plus.*plus_dummy+0*zero_dummy+I_opt_minus.*minus_dummy;
    %%%I_updated=I_opt_plus;%%%%%%%%

elseif update_spec=="gradient" % update_spec=="analytical"
    [inv_cost,inv_cost_diff]=...
        inv_cost_func(k_t,I_initial,stoch_inv_cost,theta);
    diff=-inv_cost_diff+beta_param*V_t1_diff;

    if isempty(I_min)==0
       diff=diff.*(1-(I_initial==I_min & diff<0));
    end
    if isempty(I_max)==0
        diff=diff.*(1-(I_initial==I_max & diff>0));
    end

    I_updated=I_initial+lambda_param*diff;
end

    if isempty(I_min)==0
       I_updated=max(I_updated,I_min);
    end
    if isempty(I_max)==0
       I_updated=min(I_updated,I_max);
    end

return