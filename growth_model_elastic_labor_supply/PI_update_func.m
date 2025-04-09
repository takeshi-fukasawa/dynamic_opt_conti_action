function [out,other_vars]=PI_update_func(n0,...
   Method,X0der,X0,delta,A,alpha,grid_EGM,grid,z0,z1,k0,n0_temp,c0,k1,gam,...
   nu,B,beta,n_nodes,weight_nodes,D,kdamp,n_grid,acceleration_spec)

    % Written by Takeshi Fukasawa in June 2024, based on the code of Maliar and Maliar (2013)
   

    global OPI_param geval_total feval_V_total
    global V

%==================================================================
% Method 5. Policy iteration (PI) directly updating action
%==================================================================


     c0=c0_analytical_func(n0,k0,z0,alpha,nu,gam,A,B);
     k1=k1_analytical_func(k0,n0,c0,z0,delta,A,alpha);% Compute next-period capital using budget 
                     % constraint (2) in MM(2013)

   spec_V_iter=[];
   spec_V_iter.TOL=1e-8;
   spec_V_iter.ITER_MAX=OPI_param; 
   if acceleration_spec==0
      spec_V_iter.update_spec=0;
   end

   %V=zeros(n_grid,1);
   %%%vf_coef=X0\V;%%%
   vf_coef=[];%%%%
   
   %spec_V_iter.update_spec=0;
   spec_V_iter.spectral_coef_0=1e-7;
   spec_V_iter.norm_spec=10;%% unit free
   %spec_V_iter.spectral_coef_0=spectral_coef0_param;

   [out,other_vars,iter_info_V_iter]=spectral_func(...
       @VF_Bellman_L_update_func,spec_V_iter,{V},...
    X0,n0,c0,k1,z1,gam,nu,B,beta,n_nodes,weight_nodes,vf_coef,D,kdamp);

    iter_info_V_iter.n_iter

   V_new=out{1};
   vf_coef=X0\V_new;

    n0_new=zeros(n_grid,1);
     for j=1:n_grid  % Solve for labor using eq. (18) in MM (2013)
          [n0_new(j,1),exitflag,n_iter,geval]=csolve('FOC_L_VFI',n0(j,1),[],0.000001,5,...
              k0(j,1),z0(j,1),A,alpha,gam,nu,B,beta,delta,...
              z1(j,:),n_nodes,weight_nodes,vf_coef,D);
          geval_total=geval_total+geval;   
     end

   feval_V_total=feval_V_total+iter_info_V_iter.feval*n_grid;
     
       out={n0_new};


    other_vars.k1=k1;
    other_vars.c0=c0;
    other_vars.n0=n0_new;
    other_vars.k0=k0;

return
