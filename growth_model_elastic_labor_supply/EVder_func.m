function EVder=EVder_func(k1,z1,n_nodes,weight_nodes,vf,D)
    
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
return
