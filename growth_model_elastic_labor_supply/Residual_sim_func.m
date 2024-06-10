%%% Based on the code of Maliar and Maliar (2013)
%%% Modified by Takeshi Fukasawa in June 2024

function [Mean_Residuals,Max_Residuals]  = Residuals_sim_func(k,z,A,alpha,gam,delta,sigma,rho,beta,c_coef,n_coef,D)

T     = size(k,1);   % Simulation length

% Gauss Hermite quadrature integration in the accuracy evaluation procedure:
% the number of integration nodes, their values and weights
%--------------------------------------------------------------------------
Qn = 10;             % Number of unidimensional integration nodes 
nshocks = 1;         % There is one stochastic shock
vcv = sigma^2;       % Variance covariance matrix
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(Qn,nshocks,vcv);
                     % Gauss Hermite quadrature nodes and weights;
                     % with one shock, n_nodes = Qn
                     
 for t=1:T;   
     
     k0 = k(t,1);         % Capital of period t
     z0 = z(t,1);         % Productivity level of period t 
             
     % Current period quantities in (k0,a0)
     % ------------------------------------      
     c0=Polynomial_2d([k0 z0],D)*c_coef; 
     n0=Polynomial_2d([k0 z0],D)*n_coef;

     k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);
     
     % Future period quantities in n_nodes integration nodes (k1,a1)
     % --------------------------------------------------------     
     z1 = z0.^rho.*exp(epsi_nodes); % Productivity in n_nodes nodes
     k1_dupl = ones(n_nodes,1)*k1;  % k1 is the same in n_nodes nodes       

     c1(1:n_nodes,1)=Polynomial_2d([k1_dupl z1(:,1)],D)*c_coef;
     n1(1:n_nodes,1)=Polynomial_2d([k1_dupl z1(:,1)],D)*n_coef;  
     
     % Unit-free residuals in the Euler equation
     % -----------------------------------------
     EE_LHS=c0.^(-gam);
     EE_RHS=weight_nodes'*(beta*c1(1:n_nodes,1).^(-gam).*(1-delta+alpha*A*z1(1:n_nodes,1).*k1.^(alpha-1).*n1(1:n_nodes,1).^(1-alpha)));
     Residuals(t,1)=EE_RHS./EE_LHS-1;

 end

 % Output of the accuracy evaluation test
 % --------------------------------------
 Mean_Residuals = log10(mean(mean(abs(Residuals(1:end))))); % Mean residuals 
 Max_Residuals = log10(max(max(abs(Residuals(1:end)))));    % Max  residuals
    
 
    