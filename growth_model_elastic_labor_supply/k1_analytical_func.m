function k1=k1_analytical_func(k0,n0,c0,delta,A,alpha)
    k1=(1-delta)*k0+A.*(k0.^alpha).*(n0.^(1-alpha))-c0;
end
