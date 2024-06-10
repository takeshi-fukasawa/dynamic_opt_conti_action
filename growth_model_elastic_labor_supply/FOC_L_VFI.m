function [diff_FOC_L,RHS,LHS] = FOC_L_VFI(n0,k0,z0,A,alpha,gam,nu,B,beta,delta,...
        z1,n_nodes,weight_nodes,vf_coef,D)
    
    c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
    k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);% Compute next-period capital using budget 
                     % constraint (2) in MM(2013)
    EVder=EVder_func(k1,z1,n_nodes,weight_nodes,vf_coef,D);

    LHS=-B*(1-n0).^(-nu);
    RHS=(-A).*(1-alpha).*z0.*(k0.^alpha).*(n0.^-alpha).*beta.*EVder;
    diff_FOC_L=LHS-RHS;
end

