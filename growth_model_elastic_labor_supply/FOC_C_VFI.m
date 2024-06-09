function [diff_FOC_C] = FOC_C_VFI(n0,c0,k1,k0,z0,A,alpha,gam,delta,nu,B,beta,z1,n_nodes,weight_nodes,vf,D,EVder)

EVder = 0;                  % Initialize the expected derivative of value 
                            % function

for j = 1:n_nodes           % Sum up the expected derivative of value 
                            % function across the integration nodes
    X1der = Polynomial_deriv_2d([k1 z1(:,j)],D);  
                            % Construct the derivatives of basis functions 
                            % of polynomial approximating value function at
                            % integration node j
    EVder = EVder + X1der*vf*weight_nodes(j); 
                            % Add up the weighted derivative of value 
                            % function across the nodes
end
LHS=-B*(1-n0).^(-nu);
RHS=(-A).*(1-alpha).*z0.*(k0.^alpha).*(n0.^-alpha).*beta.*EVder;
diff_FOC=LHS-RHS;

ratio=LHS./RHS;

end

