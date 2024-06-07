% Residuals_L.m is a routine for computing the residuals in the neoclassical
% growth model with elastic labor supply, considered in the 
% article "Envelope Condition Method versus Endogenous Grid Method for Solving 
% Dynamic Programming Problems" by Lilia Maliar and Serguei Maliar,  Economics 
% Letters (2013), 120, 262-266 (henceforth, MM, 2013).   
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Inputs:    "k" is current capital;
%            "z" is current productivity;
%            "A" is the normalizing constant in production;
%            "alpha" is the share of capital in production;
%            "gam, nu, B" are the utility-function parameters;
%            "delta" is the depreciation rate;
%            "sigma" is the standard deviation of the productivity shock;
%            "rho" is the autocorrelation coefficient in the process for
%                shock;
%            "beta" is the discount factor; 
%            "vf" is the vector of polynomial coefficients in the 
%               approximated value function;  
%            "D" is the degree of polynomial approximation
%
% Output:    "Mean_Residuals" and "Max_Residuals" are, respectively, the 
%            mean and maximum absolute Euler equation residuals (in log10)
% -------------------------------------------------------------------------
% Copyright © 2011-2016 by Lilia Maliar and Serguei Maliar. All rights 
% reserved. The code may be used, modified and redistributed under the 
% terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [Mean_Residuals,Max_Residuals]  = Residuals_L(k,z,A,alpha,lss,gam,nu,B,delta,sigma,rho,beta,vf,D)

T     = size(k,1);   % Simulation length

% Gauss Hermite quadrature integration in the accuracy evaluation procedure:
% the number of integration nodes, their values and weights
%--------------------------------------------------------------------------
Qn = 10;             % Number of unidimensional integration nodes 
nshocks = 1;         % There is one stochastic shock
vcv = sigma^2;       % Variance covariance matrix
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(Qn,nshocks,vcv);
                     % Gauss Hermite quadrature nodes and weights;
                     % With one shock, n_nodes = Qn
                     
 for t=1:T;   
     
     k0 = k(t,1);         % Capital of period t
     z0 = z(t,1);         % Productivity level of period t 
             
     % Current period quantities in (k0,a0)
     % ------------------------------------     
     Vder0= Polynomial_deriv_2d([k0 z0],D)*vf;  
                                  % Construct the derivative of value function
     n0 = csolve('Labor_ECM',1-lss,[],0.000001,5,nu,alpha,delta,B,A,k0,z0,Vder0);
                                   % Solve for labor using eq. (18) in MM (2013) 
     c0 = (1./(A*z0.*k0.^alpha.*n0.^-alpha.*(1-alpha)).*B.*(1-n0).^-nu).^(-1/gam);
                                   % Compute consumption using eq. (4) in MM (2013)
     k1 = (1-delta)*k0+A*z0.*k0.^alpha.*n0.^(1-alpha)-c0;                                                       
                                   % Compute next-period capital using budget 
                                   % constraint (2) in MM(2013)
   
     % Future period quantities in n_nodes integration nodes (k1,a1)
     % -------------------------------------------------------------     
     z1 = z0.^rho.*exp(epsi_nodes); % Productivity in n_nodes nodes
     k1_dupl = ones(n_nodes,1)*k1;  % k1 is the same in n_nodes nodes       
     Vder1(1:n_nodes,1) =  Polynomial_deriv_2d([k1_dupl z1(:,1)],D)*vf;   
                                    % Construct the next period derivative
                                    % of value function in the integration
                                    % nodes
     for j=1:n_nodes           
         n1(j,1) = csolve('Labor_ECM',n0,[],0.000001,5,nu,alpha,delta,B,A,k1_dupl(j,1),z1(j,1),Vder1(j,1));   
                                       % Solve for labor in all future nodes
                                       % using eq. (18) in MM (2013) 
         c1(j,1) = (1./(A*z1(j,1).*k1_dupl(j,1).^alpha.*n1(j,1).^-alpha.*(1-alpha)).*B.*(1-n1(j,1)).^-nu).^(-1/gam);
                                       % Compute consumption in all future 
                                       % nodes using eq. (4) in MM (2013)
     end  
     % Unit-free residuals in the Euler equation
     % -----------------------------------------
     Residuals(t,1) = weight_nodes'*(beta*c1(1:n_nodes,1).^(-gam)./c0.^(-gam).*(1-delta+alpha*A*z1(1:n_nodes,1).*k1.^(alpha-1).*n1(1:n_nodes,1).^(1-alpha)))-1;

 end

 % Output of the accuracy evaluation test
 % --------------------------------------
 Mean_Residuals = log10(mean(mean(abs(Residuals(1:end))))); % Mean residuals 
 Max_Residuals = log10(max(max(abs(Residuals(1:end)))));    % Maximum residuals
    

    