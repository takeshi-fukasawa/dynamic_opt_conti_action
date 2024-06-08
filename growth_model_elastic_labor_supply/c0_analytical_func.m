function c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B)
    denom=B*(1-n0).^(-nu);
    numer=A.*z0.*(1-alpha).*(k0.^alpha).*(n0.^(-alpha));
    c0=(numer./denom).^(1/gam);
end

