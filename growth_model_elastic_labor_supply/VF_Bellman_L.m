% VF_Bellman_L.m is a routine for computing value function on the left side
% of the Bellman equation in the neoclassical growth model with elastic labor 
% supply, considered in the article "Envelope Condition Method versus 
% Endogenous Grid Method for Solving Dynamic Programming Problems" by Lilia 
% Maliar and Serguei Maliar,  Economics Letters (2016), 120, 262-266 (hence-
% forth, MM, 2013).  
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Inputs:    "c0" is current consumption;
%            "n0" is current labor;
%            "k1" is next-period capital;
%            "z1" is future productivity;
%            "gam, nu, B" are the utility-function parameters;
%            "beta" is the discount factor; 
%            "n_nodes" is the integration nodes;
%            "weight_nodes" is the integration weights;
%            "vf_coef" is the vector of polynomial coefficients in the 
%               approximated value function;  
%            "D" is the degree of polynomial approximation
%
% Output:    "V_new" is value function on the left side of the Bellman
%            equation
% -------------------------------------------------------------------------
% Copyright © 2011-2016 by Lilia Maliar and Serguei Maliar. All rights 
% reserved. The code may be used, modified and redistributed under the 
% terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D)


for j = 1:n_nodes                         % For each integration node...
    X1 = Polynomial_2d([k1 z1(:,j)],D);   % Construct polynomial   
    EV(:,j) = X1*vf_coef;                 % Evaluate value function
end
if gam==1; uc0=log(c0); else uc0=(c0.^(1-gam)-1)/(1-gam); end
            % if gam=1, the consumption utility subfunction is
            % logarithmic; otherwise, it is power
if nu==1; un0=log(1-n0); else un0=((1-n0).^(1-nu)-1)/(1-nu); end
            % if nu=1, the leisure utility subfunction is
            % logarithmic; otherwise, it is power
V_new = uc0+B*un0+beta*EV*weight_nodes; % Bellman equation        
