function [out,other_vars]=update_func_L(input_cell,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,vf_coef,D,kdamp,n_grid,opts,spectral_spec)
  
  % Based on Maliar and Maliar (2013)
  % Modified by Takeshi Fukasawa in June 2024


if Method==1
    V=input_cell;
elseif Method==3
   Vder0=input_cell;
end



%==================================================================
% Method 1. Envelope condition method iterating on value function (ECM-VF)
%==================================================================
if Method==1        
     vf_coef=X0\V;
     Vder0 = X0der*vf_coef;   % Compute the derivative of value function
               
     for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
          n0(j,1)=csolve('Labor_ECM',n0(j,1),[],0.000001,5,nu,alpha,delta,B,A,k0(j,1),z0(j,1),Vder0(j,1));   
     end
                            
     c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
     k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);% Compute next-period capital using budget 
                     % constraint (2) in MM(2013)

    [V_new] = VF_Bellman_L(n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D);
                                 % Recompute value function using 
                                 % the Bellman equation

    V_new=kdamp*V_new+(1-kdamp)*V; % Update V using damping
  
elseif Method==3
%==================================================================
% Method 3. Envelope condition method iterating on derivative of value 
%           function (ECM-DVF)
%==================================================================
                      
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
            Vder0_new = kdamp*Vder0_new + (1-kdamp)*Vder0;   
                                     % Update V_new using damping
   end

            
        
%%%%%%%%
if Method==1
    out={V_new};
elseif Method==3
    out={Vder0_new};
end

other_vars.k1=k1;
other_vars.c0=c0;
other_vars.n0=n0;
other_vars.k0=k0;

return
