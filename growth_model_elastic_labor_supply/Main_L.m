% This MATLAB software solves a neoclassical stochastic growth model with
% elastic labor supply using four alternative solution methods and compares 
% their accuracy and cost. The methods considered are: (1) the envelope 
% condition method iterating on value function (ECM-VF); (2) the endogenous 
% grid method iterating on value function (EGM-VF); (3) the envelope condition 
% method iterating on derivative of value function (ECM-DVF), endogenous grid 
% method iterating on derivative of value function (EGM-DVF). For a description
% of these solution methods, see the article "Envelope Condition Method 
% versus Endogenous Grid Method for Solving Dynamic Programming Problems" by 
% Lilia Maliar and Serguei Maliar,  Economics Letters (2013), 120, 262-266 
% (henceforth, MM, 2013).  
%
% First version: July 21, 2011 
% This version:  November 2, 2016
% ------------------------------------------------------------------------
% The software uses the following files: 
% ------------------------------------------------------------------------
% 1. "Main_L.m"             solves the model and computes the residuals by 
%                           using four alternative solution methods
% 2. "VF_Bellman_L.m"       evaluates value function on the left side of
%                           the Bellman equation
% 3. "Simulation_L.m"       simulates a time-series solution given the 
%                           constructed value function
% 4. "Residuals_L.m         computes the residuals in the equilibrium conditions
%                           on a given set of points in the state space for
%                           the given numerical solution
% 5. "Polynomial_2d.m"      constructs the sets of basis functions for ordinary
%                           polynomials of the degrees from one to five, for 
%                           the model with 2 state variables
% 6."Polynomial_deriv_2d.m" constructs the derivatives of basis functions of 
%                           complete ordinary polynomial of the degrees from 
%                           one to five, for the model with 2 state variables
% 7. "GH_Quadrature.m"      constructs integration nodes and weights for the 
%                           Gauss-Hermite rules with the number of nodes in
%                           each dimension ranging from one to ten; borrowed 
%                           from Judd, Maliar and Maliar (QE, 2011)
% 8. "epsi_test.mat"        contains a fixed sequence of random numbers
% 9. "csolve.m"             a solver for a nonlinear equation (works faster
%                           than "fsolve"); written by Christopher Sims 
% -------------------------------------------------------------------------
% Copyright © 2011-2016 by Lilia Maliar and Serguei Maliar. All rights 
% reserved. The code may be used, modified and redistributed under the 
% terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

clc;
clear all;

Method = 1;   % Choose a solution method: "1", "2", "3", "4"

% "1" - envelope condition method iterating on value function (ECM-VF)
% "2" - endogenous grid method iterating on value function (EGM-VF)
% "3" - envelope condition method iterating on derivative of value 
%       function (ECM-DVF)
% "4" - endogenous grid method iterating on derivative of value 
%       function (EGM-DVF)

% 1. Model's parameters
% ---------------------
gam     = 5;                         % Utility-function parameter on consumption
nu      = 5;                         % Utility-function parameter on leisure
alpha   = 1/3;                        % Capital share in production
lss     = 2/3;                        % Steady-state share of leisure in 
                                      % total time endowment 
pkss    = 10;                         % Capital to output ratio
pcss    = 3/4;                        % Consumption to output ratio
delta   = (1-pcss)/pkss;              % Depreciation rate           
beta    = 1/(1-delta+alpha/pkss);     % Discount factor
A       = (1/beta-1+delta)/alpha/(1-lss)^(1-alpha);
                                      % Technology level that normalizes
                                      % steady state capital to one  
B       = A*(1-alpha)*(pcss/pkss)^(-gam)*lss^nu*(1-lss)^(-alpha); 
                                      % Utility-function parameter                                                                                                                                                           
rho     = 0.95;                       % Persistence of the log of the 
                                      % productivity level
sigma   = 0.01;                       % Standard deviation of shocks to the 
                                      % log of the productivity level
                                     
% 2. Steady state
% ---------------
kss = ((1/beta-1+delta)/alpha/A)^(1/(alpha-1))*(1-lss); 
                                      % Steady state capital is one 
yss = kss/pkss;                       % Steady state output                         
css = pcss*yss;                       % Steady state consumption                         

% 3. Construct a tensor product grid for computing a solution
%------------------------------------------------------------
% Unidimensional grid for capital
kmin    = 0.9;                        % Minimum capital on the grid
kmax    = 1.1;                        % Maximum capital on the grid
n_gridk = 10;                         % Number of grid point for capital 
gridk   = linspace(kmin,kmax,n_gridk)'; 
                                      % Unidimensional grid for capital of 
                                      % n_grid points

% Unidimensional grid for productivity                        
amin    = 0.9;                        % Minimum productivity on the grid
amax    = 1.1;                        % Maximum productivity on the grid
n_grida = 10;                         % Number of grid point for productivity 
grida   = linspace(amin,amax,n_grida)'; 
                                      % Unidimensional grid for productivity
% Two-dimensional tensor product grid 
grid = []; 
for i = 1:n_grida;
    grid = [grid;[gridk ones(n_gridk,1)*grida(i,1)]]; 
end

grid_EGM = grid;                     % Grid for the EGM method  

n_grid   = n_grida*n_gridk;  
                                     % Number of points in the two-dimensional 
                                     % grid 

k0 = grid(:,1);                      % Grid points for capital in the tensor 
                                     % product grid
z0 = grid(:,2);                      % Grid points for productivity in the 
                                     % tensor product grid
    
% 4. Gauss Hermite quadrature
%----------------------------
Qn      = 5;         % Number of integration nodes in Gauss Hermite quadrature
nshocks = 1;         % Number of stochastic shocks is 1
vcv     = sigma^2;   % Variance covariance matrix
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(Qn,nshocks,vcv);
                     % Gauss Hermite quadrature: number of nodes, values of 
                     % nodes and their weights, respectively
z1 = z0.^rho*exp(epsi_nodes');  
                     % Compute future shocks in all grid points and all 
                     % integration nodes; n_grid-by-n_nodes


% 5. Constructing initial guess for value function on the grid 
% ------------------------------------------------------------
D  = 2;              % Initial guess for value function is constructed using
                     % ordinary polynomial of degree 2
                     
X0   = Polynomial_2d(grid,D);
                     % Construct the matrix of explanatory variables X0 
                     % on the grid; the columns of X0 are given by basis 
                     % functions of polynomial of degree "D"        
                     
kdamp     = 1;       % Damping parameter for (fixed-point) iteration on 
                     % value function; a default value "1" means no damping
                     % and a smaller value in the interval (0,1) means that 
                     % value function is updated only partially from one
                     % iteration to another to enhance convergence
                     
vf_coef = zeros(6,1);% Initial guess for the coefficients of value function 

V = X0*vf_coef;      % Initialize value function on the grid; this grid 
                     % function will be used to check the convergence 

difv       = 1e+10;  % Initially, set the difference between value functions
                     % on two subsequent iterations to exceed the convergence
                     % criterion

n0   =  z0.*(1-lss); % Initial guess for labor
c0   =  A*z0.*k0.^alpha.*n0.^(1-alpha)*(css/yss);        
                     % Initial guess for consumption 
k1   =  k0*(1-delta)+A*z0.*k0.^alpha.*n0.^(1-alpha)-c0;  
                     % Initial guess for capital
                     
                     % Our inital guess on the labor function is that labor 
                     % changes in proportion to the productivity level; 
                     % Our initial guess on the consumption function is that 
                     % a constant fraction css/yss of the period output 
                     % goes to consumption and the rest goes to investment,
                     % where css/yss is calculated in the steady state
                                

while difv > 1e-6                          % Unit-free convergence criterion
    
        [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);

        vf_coef_new = X0\V_new;            % New vector of coefficients 
        
        vf_coef = kdamp*vf_coef_new + (1-kdamp)*vf_coef;   
                                           % Update the coefficients
        difv = max(abs(1-V_new./V));       % Compute the difference between
                                           % new and old value functions

        V = V_new;                         % Store new value function

end

% 5. Main iterative cycle: constructing polynomial approximations of degrees
% from 2 to 5 for value function 
% --------------------------------------------------------------------------

for D = 2:5;                            % For polynomial degrees from 2 to 5....   
        
    tic;
    Degree(D) = D;                      % Degree of polynomial approximation
    X0 = Polynomial_2d(grid,D);         % Ordinary polynomial of degree D 
    X0der = Polynomial_deriv_2d(grid,D);% Derivative of the polynomial
    vf_coef = X0\V;                     % Initial guess for value function
    K_coef  = X0\k1;                    % Initial guess for policy function
    difk =   inf;                       % Initial criterion of convergence 
    k_old = inf(size(grid,1),1);        % Initialize capital choices on the 
                                        % grid for checking convergence     
    
    while difk > 1e-9;                  % Convergence criterion

        %==================================================================
        % Method 1. Envelope condition method iterating on value function (ECM-VF)
        %==================================================================
        if Method==1        
       
            Vder0 = X0der*vf_coef;   % Compute the derivative of value function
               
            for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
                n0(j,1)=csolve('Labor_ECM',n0(j,1),[],0.000001,5,nu,alpha,delta,B,A,k0(j,1),z0(j,1),Vder0(j,1));   
            end
                            
            c0=(1./(A*z0.*k0.^alpha.*n0.^-alpha.*(1-alpha)).*B.*(1-n0).^-nu).^(-1/gam);
                            % Compute consumption using eq. (4) in MM (2013)

            k1 = (1-delta)*k0+A*z0.*k0.^alpha.*n0.^(1-alpha)-c0; 
                            % Compute next-period capital using budget 
                            % constraint (2) in MM(2013)

            [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                        % Recompute value function using 
                                        % the Bellman equation
            vf_coef_new = X0\V_new;     % Compute the new coefficients for 
                                        % value function       
            vf_coef = kdamp*vf_coef_new + (1-kdamp)*vf_coef;   
                                        % Update the coefficients using 
                                        % damping
        %==================================================================
        % Method 2. Endogenous grid method iterating on value function (EGM-VF)
        %==================================================================                                   
        elseif Method==2;         
            
            k1 = grid_EGM(:,1); % Grid points for next-period capital 
                                % (fixing endogenous grid)
            for j = 1:n_nodes                   
                Xder1_EGM = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1_EGM(:,j) = Xder1_EGM*vf_coef; 
                                % Compute derivative of value function in
                                % the integration nodes    
            end
            
            Wder1 = Vder1_EGM*weight_nodes;   
                                % Compute the expected derivative of next-period 
                                % value function
            for j=1:n_grid      % Solve for labor using eq. (17) in MM (2013)
                n0(j,1)=csolve('Labor_EGM',n0(j,1),[],0.000001,5,nu,gam,alpha,delta,beta,B,A,k1(j,1),z0(j,1),Wder1(j,1));
            end                                                
            
            c0 = (beta*Wder1).^(-1/gam);      
                                % Compute consumption from eq. (5) in MM
                                % (2016) 
            k0 = (B*(1-n0).^-nu./(1-alpha)/A./z0/beta./Wder1).^(1/alpha).*n0;
                                % Compute current capital from eq. (4) in
                                % MM (2013)
                                              
            [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                % Recompute value function using the Bellman 
                                % equation
                     
            grid(:,1) = k0;     % Grid points for current capital 
            
            X0 = Polynomial_2d(grid,D);   
                                % Construct polynomial on the current state 
                                % variables
            
            vf_coef_new = X0\V_new;
                                % Compute new coefficients for value function       
            vf_coef = kdamp*vf_coef_new + (1-kdamp)*vf_coef;   
                                % Update the coefficients using damping
                                   
        %==================================================================
        % Method 3. Envelope condition method iterating on derivative of value 
        %           function (ECM-DVF)
        %==================================================================
        elseif Method==3        
       
            Vder0 = X0der*vf_coef;   % Compute the derivative of value function
               
            for j=1:n_grid           % Solve for labor using eq. (18) in MM (2013)
                n0(j,1)=csolve('Labor_ECM',n0(j,1),[],0.000001,5,nu,alpha,delta,B,A,k0(j,1),z0(j,1),Vder0(j,1));   
            end
                            
            c0=(1./(A*z0.*k0.^alpha.*n0.^-alpha.*(1-alpha)).*B.*(1-n0).^-nu).^(-1/gam);
                                     % Compute consumption using eq. (4) in MM (2013)

            k1 = (1-delta)*k0+A*z0.*k0.^alpha.*n0.^(1-alpha)-c0; 
                                     % Compute next-period capital using budget 
                                     % constraint (2) in MM(2013)
                                        
            for j = 1:n_nodes           
                X1der = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1(:,j) = X1der*vf_coef; 
                                     % Compute the derivative of value
                                     % function in the integration nodes
            end
            
            Vder0_new = (1-delta+A*alpha*z0.*k0.^(alpha-1).*n0.^(1-alpha)).*(beta*Vder1*weight_nodes);
                                     % Recompute the derivative of value 
                                     % function in the grid points
            warning('off')           % Some polynomial terms are zero for
                                     % the derivative, and the system is 
                                     % underdetermined. The least-squares 
                                     % problem is still correctly 
                                     % processed by the truncated QR method
                                     % but the system produces warning
            vf_coef_new = X0der\Vder0_new;
                                     % Compute new coefficients for value 
                                     % function       
            vf_coef = kdamp*vf_coef_new + (1-kdamp)*vf_coef;   
                                     % Update the coefficients using damping
                          
        %==================================================================
        %  Method 4. Endogenous grid method iterating on derivative of value 
        %            function (EGM-DVF)
        %==================================================================
        elseif Method==4;         
            
            k1 = grid_EGM(:,1);      % Grid points for next-period capital 
                                     % (fixing endogenous grid)
            for j = 1:n_nodes              
                Xder1_EGM = Polynomial_deriv_2d([k1 z1(:,j)],D);
                Vder1_EGM(:,j) = Xder1_EGM*vf_coef; 
                                     % Compute the derivative of value function 
                                     % in the integration nodes    
            end
            
            Wder1 = Vder1_EGM*weight_nodes;   
                                     % Compute expected derivative  of
                                     % next-period value function
            for j=1:n_grid           % Solve for labor using eq. (17) in MM 
                                     % (2013)
                n0(j,1)=csolve('Labor_EGM',n0(j,1),[],0.000001,5,nu,gam,alpha,delta,beta,B,A,k1(j,1),z0(j,1),Wder1(j,1));
            end                                                
            
            c0 = (beta*Wder1).^(-1/gam);   
                                     % Compute consumption from eq. (5) in
                                     % MM (2013)
            k0 = (B*(1-n0).^-nu./(1-alpha)/A./z0/beta./Wder1).^(1/alpha).*n0;
                                     % Compute current capital from eq. (4)
                                     % in MM (2013)
                                              
            Vder0_new = (1-delta+A*alpha*z0.*k0.^(alpha-1).*n0.^(1-alpha)).*(beta*Wder1);
                                     % Recompute the derivative of value
                                     % function in grid points
            warning('off')           % Some polynomial terms are zero for
                                     % the derivative and system is 
                                     % underdetermined. The least-squares 
                                     % problem is still correctly 
                                     % processed by the truncated QR method
                                     % but the system produces warning
                   
            grid(:,1) = k0;          % Grid points for current capital 
             
            X0der = Polynomial_deriv_2d(grid,D); 
                                     % Construct polynomial on the current  
                                     % state variables
            vf_coef_new = X0der\Vder0_new;
                                     % Compute new coefficients for value function       
            vf_coef = kdamp*vf_coef_new + (1-kdamp)*vf_coef;   
                                     % Update the coefficients using damping
        %==================================================================
                          
        end               
              
        % Checking convergence 
        % --------------------
        if (Method==2)||(Method==4); 
            difk = max(abs(1-k0./k_old))     
                                    % For EGM, we check convergence of 
                                    % current capital
            k_old = k0;        
        else         
            difk = max(abs(1-k1./k_old))
                                   % For other methods, we check convergence  
                                   % of next period capital
            k_old = k1;
        end
    end
    
    % After the solution is computed by any method, we construct value 
    % function for the constructed policy rules
    
    difv = 1e10;                  % Initially, convergence criterion is not
                                  % satisfied
    X0 = Polynomial_2d(grid,D);   % Construct polynomial on the current 
                                  % state variables
                                   
    while difv > 1e-10            % Unit-free convergence criterion
        
        [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
        vf_coef_new = X0\V_new;   % New vector of coefficients 
        vf_coef = kdamp*vf_coef_new + (1-kdamp)*vf_coef;   
                                       % Update the coefficients
        difv = max(abs(1-V_new./V));   % Compute the difference between
                                       % new and old value functions
        V = V_new;                     % Store new value function
    end  

    V = V_new;                         % Update value function to be 
                                       % used as an initial guess for a 
                                       % higher degree polynomial
    CPU(D) = toc;                      % Store running time
    VK(1:1+D+D*(D+1)/2,D) = vf_coef;   % Store the solution coefficients 
       
end

% Evaluate residuals in the model's equations on a simulated path
% ---------------------------------------------------------------
fprintf(1,'Method = %i:\n',Method);
fprintf(1,'ACCURACY EVALUATION AND RUNNING TIME:\n\n');
for D = 2:5 % For polynomial degrees from 2 to 5... 
    [k,z] = Simulation_L(A,alpha,lss,gam,nu,B,delta,sigma,rho,VK(1:1+D+D*(D+1)/2,D),D);
            % Simulate the solution under a sequence of 10,200 shocks stored
            % "epsi_test.mat" and discard the first 200 entries to eliminate
            % the effect of initial conditions
    [Mean_Residuals(D),Max_Residuals(D)] = Residuals_L(k,z,A,alpha,lss,gam,nu,B,delta,sigma,rho,beta,VK(1:1+D+D*(D+1)/2,D),D);
            % Compute residuals in the model's equations in each point of the 
            % simulated path and compute the mean and maximum residuals in
            % the model's equations
    fprintf(1,'Polynomial of degree = %i:\nRunning time = %.2f, Mean residuals = %.2f, Max residuals = %.2f\n\n',Degree(D),CPU(D),Mean_Residuals(D),Max_Residuals(D));
            % Display the results
end
