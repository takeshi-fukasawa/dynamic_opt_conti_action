function [diff_FOC_C,RHS,LHS] = FOC_C_VFI(EVder,n0,c0,k0,z0,A,alpha,gam,nu,B,beta)

    LHS=c0.^(-gam);
    RHS=beta*EVder;
    diff_FOC_C=LHS-RHS;

end

