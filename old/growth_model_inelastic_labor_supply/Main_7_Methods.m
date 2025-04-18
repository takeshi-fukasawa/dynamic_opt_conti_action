function [out,other_output]=growth_main(Method,spectral_spec,D)
% This MATLAB software solves a neoclassical stochastic growth model using
% seven alternative solution methods and compares their accuracy and cost: 
% envelope condition method iterating on value function (ECM-VF), conventional 
% value function iteration (VFI), endogenous grid method (EGM), policy 
% iteration using envelope condition (ECM-PI), conventional policy iteration 
% (ECM-PI), envelope condition method iteration on the derivative of value 
% function (ECM-DVF) and conventional Euler equation method (EE). For a 
% description of these solution methods, see the article "Envelope Condition 
% Method with an Application to Default Risk Models" by Cristina Arellano, 
% Lilia Maliar, Serguei Maliar and Viktor Tsyrennikov, Journal of Economic 
% Dynamics and Control (2016), 69, 436-459 (henceforth, AMMT, 2016).  
%
% First version: December 12, 2012 
% This version:  October 31, 2016
% Modified by Takeshi Fukasawa in May 2024

% ------------------------------------------------------------------------
% The software uses the following files: 
% ------------------------------------------------------------------------
% 1. "Main_7.m"             solves the model and computes the residuals by 
%                           using seven alternative solution methods
% 2. "VF_Bellman.m"         evaluates value function on the left side of
%                           the Bellman equation
% 3. "FOC_VFI.m"            computes new capital that satisfies the FOC 
% 4. "EGM_BC.m"             computes current capital that satisfies budget 
%                           constraint given next-period capital 
% 5. "Simulation_ECM.m"     simulates a time-series solution given the 
%                           constructed value function
% 6. "Residuals_ECM.m       computes the residuals in equilibrium conditions
%                           on a given set of points in the state space for
%                           a given numerical solution
% 7. "Polynomial_2d.m"      constructs the sets of basis functions for ordinary
%                           polynomials of the degrees from one to five, for 
%                           the model with 2 state variables
% 9."Polynomial_deriv_2d.m" constructs the derivatives of basis functions of 
%                           complete ordinary polynomial of the degrees from 
%                           one to five, for the model with 2 state variables
% 9. "GH_Quadrature.m"      constructs integration nodes and weights for the 
%                           Gauss-Hermite rules with the number of nodes in
%                           each dimension ranging from one to ten; borrowed 
%                           from Judd, Maliar and Maliar (QE, 2011)
% 10. "epsi_test.mat"       contains a fixed sequence of random numbers
% -------------------------------------------------------------------------
% Copyright � 2012-2016 by Lilia Maliar and Serguei Maliar. All rights 
% reserved. The code may be used, modified and redistributed under the 
% terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

fprintf('\n\n\n\n\nBeginning execution with method %i\n', Method)

% Method = 1;   % Choose a solution method: "1", "2", "3", "4", "5", "6", "7"

% "1" - envelope condition method iterating on value function (ECM-VF)
% "2" - conventional value function interation (VFI)
% "3" - endogenous grid method (EGM)
% "4" - policy function iteration via envelope condition (ECM-PI) 
% "5" - conventional policy function iteration via FOC (PI)
% "6" - envelope condition method iterating on derivative of value function (ECM-DVF)
% "7" - conventional Euler equation method (EE)

global iter_info iter_info0 V k1
global spectral_coef0_param lambda_param
global common_spectral_coef_spec

D_init=D;
D_min=D;
D_max=D;

% 1. Model's parameters
% ---------------------
gam     = 2.0;       % Utility-function parameter: the code is stable for 
                     % the range of [0.08, 8]
alpha   = 0.36;      % Capital share in output
beta    = 0.99;      % Discount factor
delta   = 0.02;      % Depreciation rate 
rho     = 0.95;      % Persistence of the log of the productivity level
sigma   = 0.01;      % Standard deviation of shocks to the log of the 
                     % productivity level
                     
% 2. Technology level
% -------------------
A       = (1/beta-(1-delta))/alpha;  
                     % Normalize steady state capital to one             
css     = A-delta;   % Steady state consumption 
yss     = A;         % Steady state output


% 3. Construct a tensor product grid for computing a solution
%------------------------------------------------------------
% Unidimensional grid for capital
kmin    = 0.9;          % min capital on the grid
kmax    = 1.1;          % max capital on the grid
n_gridk = 10;           % number of grid points for capital 
gridk   = linspace(kmin,kmax,n_gridk)'; 
                        % Unidimensional grid for capital of n_grid points

% Unidimensional grid for productivity                        
zmin    = 0.9;          % min productivity on the grid
zmax    = 1.1;          % max productivity on the grid
n_gridz = 10;           % number of grid points for productivity 
gridz   = linspace(zmin,zmax,n_gridz)'; 
                        % Unidimensional grid for productivity

% Two-dimensional tensor product grid 
grid = []; 
for i = 1:n_gridz;
    grid = [grid;[gridk ones(n_gridk,1)*gridz(i,1)]]; 
end

grid_EGM = grid;        % Grid for the EGM method  

n_grid   = n_gridz*n_gridk;  
                        % Number of points in the two-dimensional grid 

k0 = grid(:,1);         % Grid points for capital in the tensor product grid
z0 = grid(:,2);         % Grid points for productivity in the tensor product grid
    
% 4. Gauss Hermite quadrature
%----------------------------
Qn      = 5;         % Number of integration nodes in Gauss Hermite quadrature
nshocks = 1;         % Number of stochastic shocks is 1
vcv     = sigma^2;   % Variance covariance matrix
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(Qn,nshocks,vcv);
                     % Gauss Hermite quadrature: number of nodes, values of 
                     % nodes and their weights, respectively
z1 = z0.^rho*exp(epsi_nodes');  
            % Compute future shocks in all grid points and all integration
            % nodes; the size is n_grid-by-n_nodes; 


% 5. Constructing initial guess for value function on the grid 
% ------------------------------------------------------------
D  = D_init;              % Initial guess for value function is constructed using
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
                     
vf_coef = zeros(size(X0,2),1);% Initial guess for the coefficients of value function 

V = X0*vf_coef;      % Initialize the value function on the grid; this
                     % grid function will be used to check the convergence 

difv       = 1e+10;  % Initially, set the difference between the value 
                     % functions on two subsequent iterations to exceed 
                     % the convergence criterion

c0   =  A*z0.*k0.^alpha*(css/yss);        % Initial guess for consumption 
k1   =  k0*(1-delta)+A*z0.*k0.^alpha-c0;  % Initial guess for capital
                     % Our initial guess for consumption function is that 
                     % a constant fraction css/yss of the period output 
                     % goes to consumption and the rest goes to investment,
                     % where css/yss is calculated in the steady state

spec_default.norm_spec=10;%% unit free
spec_default.TOL=1e-6;
spec_default.ITER_MAX=2000;
spec_default.spectral_coef_0=spectral_coef0_param;
spec_default.common_spectral_coef_spec=common_spectral_coef_spec;

spec=spec_default;
spec.ITER_MAX=2000;
 
    [output_spectral,other_vars,iter_info_V]=...
        spectral_func(@VF_Bellman_update_func,spec,{V},...
            X0,c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

    V=output_spectral{1};

    if max(abs(V))>10^10
       iter_info_V.FLAG_ERROR=1;
    end

    vf_coef = X0\V;     % Coefficients for value function

%%% Initial values %%%
V_init=V;
vf_coef_init=vf_coef;
%%%%%%%%%%

% 5. Main iterative cycle: constructing polynomial approximations of degrees 
% from 2 to 5 for value function 
% --------------------------------------------------------------------------

for D = D_min:D_max;                            % For polynomial degrees from 2 to 5...   
        
    tStart=tic;
    Degree(D) = D;                      % Degree of polynomial approximation
    X0 = Polynomial_2d(grid,D);         % Ordinary polynomial of degree D 
    X0der = Polynomial_deriv_2d(grid,D);% Derivative of the polynomial
    vf_coef = X0\V;                     % Initial guess for value function
    K_coef  = X0\k1;                    % Initial guess for policy function
    difk =   inf;                       % Initial criterion of convergence 
    k_old = inf(size(grid,1),1);        % Initialize capital choices on the 
                                        % grid for checking convergence 
    opts = optimset('Display','none','Algorithm','trust-region-dogleg','MaxFunEvals',10000,'MaxIter',1000,'TolX',1e-10,'TolFun',1e-10);            
                                        % Options for the solver
    Vder0 = X0der*vf_coef;              % Compute derivative of value function


    spec=spec_default;
    spec.TOL=1e-9;

    if spectral_spec==0
        spec.update_spec=0;
    elseif spectral_spec==2;
        spec.SQUAREM_spec=1;
    end

    if Method==1 | Method==2
        input={V};
    elseif Method==3
        input={vf_coef};
    elseif Method==4 | Method==5 | Method==7
        input={k1};
    elseif Method==6
        input={Vder0};
    elseif Method==0
        input={V,k1};
       TOL_vec=(spec.TOL)*ones(1,2);
       TOL_vec(2)=TOL_vec(2)*lambda_param;%%% TOL of action should not depend on lambda_param 
    end

    if Method<=7 & Method>=1
        fun=@update_func;
    elseif Method==0
        fun=@joint_update_func;
    end

    [output_spectral,other_vars,iter_info]=...
        spectral_func(fun,spec,input,...
            Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,k1,gam,c0,beta,...
            n_nodes,weight_nodes,vf_coef,D,kdamp,n_grid,opts,spectral_spec);
    %iter_info.feval;
    

    if Method==1 | Method==2
        V=output_spectral{1};
        vf_coef = X0\V;     % Coefficients for value function 
        k1=other_vars.k1;
    elseif Method==3
        vf_coef=output_spectral{1};
        k1=other_vars.k1;
        k0=other_vars.k0;
        
        grid(:,1) = k0;        % Grid points for current capital  
        X0 = Polynomial_2d(grid,D);   % Construct polynomial on 
                                          % current state variables

        vf_coef=X0\V; % Coefficients for value function

    elseif Method==4 | Method==5 | Method==7
        k1=output_spectral{1};
        V=other_vars.V;
        vf_coef = X0\V;     % Coefficients for value function 
    
    elseif Method==6
       Vder0=output_spectral{1};
        vf_der_coef = X0\Vder0;     % Coefficients for value function
        k1=other_vars.k1;

    elseif Method==0
        V=output_spectral{1};
        k1=output_spectral{2};
        vf_coef=X0\V;
    end


    c0=other_vars.c0;
    k0=other_vars.k0;

    
    k_coef = X0\k1;
    c0=(1-delta)*k0+A*z0.*k0.^alpha-k1;
    c_coef = X0\c0;

    if max(abs(k1))>10^10 || iter_info.feval==iter_info.ITER_MAX
       iter_info.FLAG_ERROR=1;
    end
    
    % After the solution is computed by any method, we construct the value 
    % value function for the constructed policy rules
    if 1==0
        spec=spec_default;
        spec.TOL=1e-10;
        [output_spectral,other_vars,iter_info_V]=...
            spectral_func(@VF_Bellman_update_func,spec,{V},...
            X0,c0,k1,z1,gam,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

        V=output_spectral{1};
        vf_coef = X0\V;     % Coefficients for value function
        
        if max(abs(V))>10^10
            iter_info_V.FLAG_ERROR=1;
        end
        feval_V(D)=iter_info_V.feval;
    else
        feval_V(D)=0;
    end

    CPU(D) = toc(tStart);                      % Store running time
    V_coef(1:1+D+D*(D+1)/2,D) = vf_coef;   % Store the solution coefficients (V)
    K_coef(1:1+D+D*(D+1)/2,D) = k_coef;   % Store the solution coefficients (V)
    C_coef(1:1+D+D*(D+1)/2,D) = c_coef;   % Store the solution coefficients (V)

    feval(D)=iter_info.feval;
    
    
end

% Evaluate residuals in the model's equations on a simulated path
% ---------------------------------------------------------------
fprintf(1,'Method = %i:\n',Method);
fprintf(1,'ACCURACY EVALUATION AND RUNNING TIME:\n\n');
for D = D_min:D_max % For polynomial degrees from 2 to 5... 
    if 1==0
        %%% Use V_coef for validating the accuracy of the solution
        [k,z] = Simulation_ECM(A,alpha,gam,delta,sigma,rho,V_coef(1:1+D+D*(D+1)/2,D),D);
            % Simulate the solution under a sequence of 10,200 shocks stored
            % "epsi_test.mat" and discard the first 200 entries to eliminate
            % the effect of initial conditions
        [Mean_Residuals(D),Max_Residuals(D)] = Residuals_ECM(k,z,A,alpha,gam,delta,sigma,rho,beta,V_coef(1:1+D+D*(D+1)/2,D),D);
            % Compute residuals in the model's equations in each point of the 
            % simulated path and compute the mean and maximum residuals in
            % the model's equations
    else
        %%% Use C_coef, K_coef for validating the accuracy of the solution
        [k,z] = Simulation_ECM2(A,alpha,gam,delta,sigma,rho,K_coef(1:1+D+D*(D+1)/2,D),D);
            % Simulate the solution under a sequence of 10,200 shocks stored
            % "epsi_test.mat" and discard the first 200 entries to eliminate
            % the effect of initial conditions
        [Mean_Residuals(D),Max_Residuals(D)] = Residuals_ECM2(k,z,A,alpha,gam,delta,sigma,rho,beta,C_coef(1:1+D+D*(D+1)/2,D),D);
            % Compute residuals in the model's equations in each point of the 
            % simulated path and compute the mean and maximum residuals in
            % the model's equations
    %end    

     fprintf(1,'Polynomial of degree = %i:\nRunning time = %.3f, Mean residuals = %.2f, Max residuals = %.2f\n\n',Degree(D),CPU(D),Mean_Residuals(D),Max_Residuals(D));
            % Display the results
    
    out(D-1,:)=[D,Degree(D),CPU(D),Mean_Residuals(D),Max_Residuals(D),...
        feval(D),feval_V(D),1-iter_info.FLAG_ERROR];

    other_output.iter_info=iter_info;
    other_output.iter_info_V=iter_info_V;
    
end

end  % function