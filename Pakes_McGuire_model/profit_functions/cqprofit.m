function [] = cqprofit(nfirms, descn)
    % competition
    global profit agprof csurplus share pmcmarg concent
    global egw M mc w sigma p

     i = 1;
     while i <= descn;
       progress(i);
       w = qdecode(i,nfirms+1);
        %p = newton(p,@cfunk);
        opts = optimset('Display','none','Algorithm','trust-region-dogleg','MaxFunEvals',10000,'MaxIter',1000,'TolX',1e-10,'TolFun',1e-10);
        p=fsolve(@cfunk,p,opts);
       profstar = M*p.*sigma - M*mc*sigma;
       profit(i,:) = profstar';
       csurplus(i) = M*log(1+sum(exp(log(egw)-p)));
       agprof(i,:) = profstar';
       share(i,:) = sigma';
       pmcmarg(i) = sum(p.*sigma) / mc / sum(sigma);
       concent(i) = max(sigma)/sum(sigma);
       i = i+1;
     end
end
