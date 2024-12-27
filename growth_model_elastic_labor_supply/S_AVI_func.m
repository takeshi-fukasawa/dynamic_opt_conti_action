function [V_sol,other_vars,iter_info]=S_AVI_func(Method,V_0,TOL,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid)
 
%==================================================================
% Method 6,7. Safe Accelerated Value Iteration (S-AVI; Goyal and Grand-Clement 2021)
%==================================================================

gamma_SAVI=(1-sqrt(1-beta^2))/beta;
alpha_SAVI=1/(1+beta);
lambda_SAVI=(1+beta)/2;
Method_VFI=0;
spectral_spec=0;
ITER_MAX=500;
safe_spec=0;%%%%
feval=0;

if Method==6
    safe_spec=0;
elseif Method==7
    safe_spec=1;
end

 [out,other_vars]=update_func_L(V_0,...
   Method_VFI,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,spectral_spec);

V_s=out{1};feval=feval+1;
V_s_minus_1=V_s;
T_V_0=V_s;
DIST_0_temp=max(abs(T_V_0(:)-V_0(:)));

for s=1:ITER_MAX
   h_s=V_s+gamma_SAVI*(V_s-V_s_minus_1);

    [out,other_vars]=update_func_L(h_s,...
      Method_VFI,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
      nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,spectral_spec);

   T_h_s=out{1};feval=feval+1;
   V_s_plus_1_temp=h_s-alpha_SAVI*(h_s-T_h_s);

   if safe_spec==1
           [out,other_vars]=update_func_L(V_s_plus_1_temp,...
            Method_VFI,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
            nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,spectral_spec);
    
           T_V_s_plus_1_temp=out{1};feval=feval+1;

           DIST_temp=max(abs(T_V_s_plus_1_temp(:)-V_s_plus_1_temp(:)));

        if DIST_temp<=DIST_0_temp*lambda_SAVI^(s+1)
              V_s_plus_1=V_s_plus_1_temp;
        else

             [out,other_vars]=update_func_L(h_s,...
                Method_VFI,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
                nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,spectral_spec);

              T_V_s=out{1};feval=feval+1;

              V_s_plus_1=T_V_s;
         end

    elseif safe_spec==0
       V_s_plus_1=V_s_plus_1_temp;
    end
    
    DIST=(max(abs((V_s_plus_1(:)-V_s(:))./V_0(:))));
    if DIST<TOL
        break;
    end

    V_s_minus_1=V_s;
    V_s=V_s_plus_1;

end%s=1:ITER_MAX

V_sol=V_s_plus_1;
iter_info.n_iter=s;
iter_info.feval=feval;
iter_info.ITER_MAX=ITER_MAX;

iter_info.FLAG_ERROR=0;
if s==ITER_MAX
    iter_info.FLAG_ERROR=1;
end

end
