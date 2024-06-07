% Simulation_L.m is a routine for simulating the solution to the
% neoclassical growth model with elastic labor supply, considered in the 
% article "Envelope Condition Method versus Endogenous Grid Method for Solving 
% Dynamic Programming Problems" by Lilia Maliar and Serguei Maliar,  Economics 
% Letters (2013), 120, 262-266 (henceforth, MM, 2013).    
% -------------------------------------------------------------------------
% Inputs:    "A" is the normalizing constant in production;
%            "alpha" is the share of capital in production;
%            "gam, nu, B" are the utility-function parameters;
%            "delta" is the depreciation rate;
%            "sigma" is the standard deviation of the productivity shock;
%            "rho" is the autocorrelation coefficient in the process for
%               shock;
%            "vf" is the vector of polynomial coefficients in the 
%               approximated value function;  
%            "D" is the degree of polynomial approximation
%
% Output:    "k" is current capital;
%            "z" is current level of productivity
% -------------------------------------------------------------------------
% Copyright © 2011-2016 by Lilia Maliar and Serguei Maliar. All rights 
% reserved. The code may be used, modified and redistributed under the 
% terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------


function [k,z] = Simulation_L(A,alpha,lss,gam,nu,B,delta,sigma,rho,vf,D)

load epsi_test;            % Load innovations; they are normally distributed 
                           % with zero mean and unit standard deviation of 
                           % length 10,200 
T     = size(epsi_test,1); % Simulation length

% Simulate the productivity series
% --------------------------------
epsi  = epsi_test*sigma;   % Standard deviation is sigma
z(1,1)  = 1;               % Initial productivity is 1
for t = 2:T; 
   z(t,1) = z(t-1,1).^rho.*exp(epsi(t,1)); % Simulation AR(1) process
end;

% Simulate the endogenous model's variables
% -----------------------------------------
k(1,1)  = 1;               % Initial capital is 1 
for t = 1:T 
    Vder(t,1) = Polynomial_deriv_2d([k(t,1) z(t,1)],D)*vf;
                           % Construct the derivative of value function
    n(t,1)=csolve('Labor_ECM',1-lss,[],0.000001,5,nu,alpha,delta,B,A,k(t,1),z(t,1),Vder(t,1));
                           % Solve for labor using eq. (18) in MM (2013) 
    c(t,1)=(1./(A*z(t,1).*k(t,1).^alpha.*n(t,1).^-alpha.*(1-alpha)).*B.*(1-n(t,1)).^-nu).^(-1/gam);
                           % Compute consumption using eq. (4) in MM (2013)
    k(t+1,1) = (1-delta)*k(t,1)+A*z(t,1).*k(t,1).^alpha.*n(t,1).^(1-alpha)-c(t,1);                                                       
                           % Compute next-period capital using budget 
                           % constraint (2) in MM(2013)
end   

% Discard initial observations
% ----------------------------
discard = 200;             % The number of observations to discard in order  
                           % to eliminate the effect of initial conditions
z       = z(discard+1:T,1);
k       = k(discard+1:T,1);