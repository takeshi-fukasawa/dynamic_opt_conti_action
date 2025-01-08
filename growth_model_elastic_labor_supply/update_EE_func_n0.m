function [out,other_vars]=update_EE_func_n0(var,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec)

    global analytical_EE_spec

    % Written by Takeshi Fukasawa in June 2024, based on the code of Maliar and Maliar (2013)
   
    n0=var;

    c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
    k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);
    
    n_coef=X0\n0;

    dEV_dk1=0;
    for j = 1:n_nodes    % Sum up the expected derivative of value
            % function across the integration nodes
            X1_j = Polynomial_2d([k1 z1(:,j)],D);  
            % Construct the derivatives of basis functions 
            % of polynomial approximating value function at
            % integration node j

            n1_j=X1_j*n_coef;

            numer=B*(1-n1_j).^(-nu);
            denom=(1-alpha).*A.*z1(:,j).*(k1.^alpha).*(n1_j.^(-alpha));

            du1_dc1=numer./denom;

            dEV_dk1=dEV_dk1+du1_dc1.*(1-delta+z1(:,j).*A.*alpha.*(k1.^(alpha-1)).*(n1_j.^(1-alpha))).*weight_nodes(j);
    end

    dk1_dn0=z0.*A.*(k0.^alpha).*(1-alpha).*(n0.^(-alpha));
    dEV_dn0=dk1_dn0.*dEV_dk1;

    if analytical_EE_spec==1% Analytical sol
        numer=B;
        denom=dEV_dn0*beta;
        n0_new=1-(numer./denom).^(1/nu);
    else % Numerically solve sol
        for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
          [n0_new(j,1),exitflag,n_iter,geval]=csolve('EE_resid_func',...
              n0(j,1),[],0.000001,5,...
              dEV_dn0(j),B,beta,nu);  
        end
    end

       out={n0_new};


    other_vars.k1=k1;
    other_vars.c0=c0;
    other_vars.n0=n0;
    other_vars.k0=k0;

return
