function [out,other_output]=Main_function(Method,spectral_spec,D)

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
% Modified by Takeshi Fukasawa in June 2024.

%clc;
%clear all;


% "1" - envelope condition method iterating on value function (ECM-VF)
% "2" - endogenous grid method iterating on value function (EGM-VF)

fprintf('\n\n\n\n\nBeginning execution with method %i\n', Method)

global iter_info iter_info0 V k1
global alpha0_param lambda_param
global common_spectral_coef_spec n0 c0

D_init=D;
D_min=D;
D_max=D;

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
      
%%%%%%%%%%%%%%%%%%
spec_default.norm_spec=10;%% unit free
spec_default.TOL=1e-6;
spec_default.ITER_MAX=3000;
spec_default.alpha_0=alpha0_param;
spec_default.common_spectral_coef_spec=common_spectral_coef_spec;
spec_default.DEBUG=1;%%%%%%%%%%%%%%%%%%%

spec=spec_default;
 
    [output_spectral,other_vars,iter_info_V]=...
        spectral_func(@VF_Bellman_L_update_func,spec,{V},...
            X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

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
    Vder0 = X0der*vf_coef;              % Compute derivative of value function


    spec=spec_default;
    spec.TOL=1e-6;

    if spectral_spec==0
        spec.update_spec=0;
    elseif spectral_spec==2;
        spec.SQUAREM_spec=1;
    end
    
    if Method==1 | Method==3 | Method==0
        input={V};
    elseif Method==2
        input={vf_coef};
    elseif Method==4
         input={V,k0};
        spec.common_spectral_coef_spec=1; % Important
    elseif Method==-1
        input={V,n0};
       TOL_vec=(spec.TOL)*ones(1,2);
       TOL_vec(2)=TOL_vec(2)*lambda_param;%%% TOL of action should not depend on lambda_param 
       spec.TOL=TOL_vec;
       spec.norm_spec=[10,0];% unit free norm for V
        
    elseif Method==-2
        input={V,n0,c0};
        TOL_vec=(spec.TOL)*ones(1,3);
        TOL_vec(1,2:3)=TOL_vec(1,2:3).*lambda_param;%%% TOL of action should not depend on lambda_param 
        spec.TOL=TOL_vec;
        spec.norm_spec=[10,0,0];% unit free norm for V

    elseif Method==-3
        input={V,[n0;c0]};
        TOL_vec=(spec.TOL)*ones(1,2);
        TOL_vec(2)=TOL_vec(2)*lambda_param;%%% TOL of action should not depend on lambda_param 
        spec.TOL=TOL_vec;
        spec.norm_spec=[10,0];% unit free norm for V
    end


    if Method==1 | Method==0
        fun=@update_func_L;
    elseif Method==2
        fun=@update_func_L;
    elseif Method==-1
        fun=@VF_PGI_func_n0;
    elseif Method==-2
        fun=@VF_PGI_func_n0_c0;
    elseif Method==-3
        fun=@VF_PGI_func_n0_c0_2;
    end


    [output_spectral,other_vars,iter_info]=...
        spectral_func(fun,spec,input,...
        Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0,c0,k1,gam,...
        nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,spectral_spec);
    %iter_info.feval;
    

    if Method==1 | Method==0
        V=output_spectral{1};
        vf_coef = X0\V;     % Coefficients for value function 
        k1=other_vars.k1;
    elseif Method==2
        vf_coef=output_spectral{1};
        k1=other_vars.k1;
        k0=other_vars.k0;
        
        grid(:,1) = k0;        % Grid points for current capital  
        X0 = Polynomial_2d(grid,D);   % Construct polynomial on 
                                          % current state variables
    elseif Method==-1
        V=output_spectral{1};
        n0=output_spectral{2};
        k1=other_vars.k1;
        vf_coef=X0\V;
    elseif Method==-2
        V=output_spectral{1};
        n0=output_spectral{2};
        c0=output_spectral{3};
        k1=other_vars.k1;
        vf_coef=X0\V;
    elseif Method==-3
        V=output_spectral{1};
        var2=output_spectral{2};
        n0=var2(1:n_grid,1);
        c0=var2(n_grid+1:n_grid*2,1);
        k1=other_vars.k1;
        vf_coef=X0\V;
    end

    c0=other_vars.c0;
    n0=other_vars.n0;
    k0=other_vars.k0;

    %%%%%%%%%%%%%%%
       
    vf_coef=X0\V; % Coefficients for value function
    EVder=EVder_func(k1,z1,n_nodes,weight_nodes,vf_coef,D);
    FOC_val_C=FOC_C(EVder,n0,c0,k0,z0,A,alpha,gam,nu,B,beta);    
    %%%%%%%%%%%%%%%

    k_coef = X0\k1;
    c_coef = X0\c0;
    n_coef = X0\n0;

    if max(abs(k1))>10^10 || iter_info.feval==iter_info.ITER_MAX
       iter_info.FLAG_ERROR=1;
    end
    
    % After the solution is computed by any method, we construct the value 
    % value function for the constructed policy rules
    if 1==0
        spec=spec_default;
        spec.TOL=1e-9;
        [output_spectral,other_vars,iter_info_V]=...
            spectral_func(@VF_Bellman_L_update_func,spec,{V},...
            X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

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
    N_coef(1:1+D+D*(D+1)/2,D) = n_coef;   % Store the solution coefficients (V)

    feval(D)=iter_info.feval;
    
    
end

% Evaluate residuals in the model's equations on a simulated path
% ---------------------------------------------------------------
fprintf(1,'Method = %i:\n',Method);
fprintf(1,'ACCURACY EVALUATION AND RUNNING TIME:\n\n');


for D = D_min:D_max % For polynomial degrees from 2 to 5... 
        %%% Use C_coef, K_coef for validating the accuracy of the solution
        [k,z] = simulation_func(A,alpha,gam,delta,sigma,rho,K_coef(1:1+D+D*(D+1)/2,D),D);             % Simulate the solution under a sequence of 10,200 shocks stored
            % "epsi_test.mat" and discard the first 200 entries to eliminate
            % the effect of initial conditions

        [Mean_Residuals(D),Max_Residuals(D)]  = Residual_sim_func(k,z,A,alpha,gam,delta,sigma,rho,beta,...
            C_coef(1:1+D+D*(D+1)/2,D),N_coef(1:1+D+D*(D+1)/2,D),D);
    

     fprintf(1,'Polynomial of degree = %i:\nRunning time = %.3f, Mean residuals = %.2f, Max residuals = %.2f\n\n',Degree(D),CPU(D),Mean_Residuals(D),Max_Residuals(D));
            % Display the results
    
    out(D-1,:)=[D,Degree(D),CPU(D),Mean_Residuals(D),Max_Residuals(D),...
        feval(D),feval_V(D),1-iter_info.FLAG_ERROR,CPU(D)/feval(D)];

    other_output.iter_info=iter_info;
    other_output.iter_info_V=iter_info_V;
    
end

end  % function