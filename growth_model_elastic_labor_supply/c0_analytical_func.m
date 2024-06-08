function c0=c0_analytical_func(n0,k0,alpha,nu,gam,A,B)
    numer=B*(1-n0).^(-nu);
    denom=A.*(1-alpha).*(k0.^alpha).*(n0.^(-alpha));
    c0=(numer./denom).^(1/gam);
end
