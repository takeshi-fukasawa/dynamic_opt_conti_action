% Simulation_ECM.m is a routine for simulating the solution to the
% neoclassical growth model, considered in the article  "Envelope Condition 
% Method with an Application to Default Risk Models" by Cristina Arellano, 
% Lilia Maliar, Serguei Maliar and Viktor Tsyrennikov, Journal of Economic
% Dynamics and Control (2016), 69, 436-459 (henceforth, AMMT, 2016).    
% -------------------------------------------------------------------------
% Inputs:    "A" is the normalizing constant in production;
%            "alpha" is the share of capital in production;
%            "gam" is the utility-function parameter;
%            "delta" is the depreciation rate;
%            "sigma" is the standard deviation of the productivity shock;
%            "rho" is the autocorrelation coefficient in the process for
%               shock;
%            "vf" is the vector of polynomial coefficients in the 
%               approximated value function;  
%            "D" is the degree of polynomial approximation
%
% Output:    "k" is current capital;
%            "z" is current productivity
% -------------------------------------------------------------------------
% Copyright � 2012-2016 by Lilia Maliar and Serguei Maliar. All rights 
% reserved. The code may be used, modified and redistributed under the 
% terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------


function [k,z] = Simulation_ECM2(A,alpha,gam,delta,sigma,rho,k_coef,D)
epsi_test = csvread(string('EE_SHOCKS.csv'));
% load epsi_test;            % Load innovations normally distributed with zero 
                             % mean and unit standard deviation of length 10,200 
T     = size(epsi_test,1); % Simulation length

% Simulate the productivity series
% --------------------------------
epsi  = epsi_test*sigma;   % Standard deviation is sigma
z(1,1)  = 1;               % Initial productivity is 1
for t = 2:T; 
   z(t,1) = z(t-1,1).^rho.*exp(epsi(t,1)); % Simulating AR(1) process
end;

% Simulate endogenous model's variables
% -------------------------------------
k(1,1)  = 1;                % Initial capital is 1 
for t = 1:T 
    k(t+1,1) = Polynomial_2d([k(t,1) z(t,1)],D)*k_coef;
end   

% Discard initial observations
% ----------------------------
discard = 200;             % The number of observations to discard to 
                           % eliminate the effect of initial conditions
z       = z(discard+1:T,1);
k       = k(discard+1:T,1);
