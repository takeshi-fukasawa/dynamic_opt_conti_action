function [out] =func_for_GMRES(V)

vf_coef=X0\V;

for j = 1:n_nodes                         % For each integration node...
    X1 = Polynomial_2d([k1 z1(:,j)],D);   % Construct polynomial   
    EV(:,j) = X1*vf_coef;                 % Evaluate value function
end

out=V-beta*EV*weight_nodes;        
