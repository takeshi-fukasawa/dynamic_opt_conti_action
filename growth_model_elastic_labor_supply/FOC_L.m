function [diff_FOC_L,RHS,LHS] = FOC_L(EVder,n0,c0,k0,z0,A,alpha,gam,nu,B,beta)
    
    LHS=-B*(1-n0).^(-nu);
    RHS=(-A).*(1-alpha).*z0.*(k0.^alpha).*(n0.^-alpha).*beta.*EVder;
    diff_FOC_L=LHS-RHS;
end

