function [resid] = EE_resid_func(n0,dEV_dn0,B,beta,nu)

    resid=-B.*(1-n0).^(-nu)+beta.*dEV_dn0;
end
