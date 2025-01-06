function [out,other_vars]=...
    policy_iter_func(I_t,V_t,...
    k_t,exo_t,pi_mat,exo_t1,basis_t,inv_multiply_t,basis_exo_t1,...
    x_inv,w_inv,...
    state_min,state_max,Smol_elem,mu_max,d,ind,parameters,I_min,I_max)

% I_t: n_pts*N*n_node_inv
% V_t: n_pts*N
% k_t: n_pts*N
% exo_t1: n_pts*N*n_node_inv
% basis_t: n_pts*n_coef
% inv_multiply_t: n_pts*n_coef
% basis_exo_t1: n_pts*n_coef*1
% x_inv: n_node_inv*1
% w_inv:n_node_inv*1

%% solve optimal I (modified version)

global beta_param delta_param spec update_spec
global spec_precompute diff
global geval_total veval_total
global OPI_param

[n_pts,N,n_node_inv]=size(I_t);
n_node_inv=size(w_inv,1);

theta=parameters(1:4);
%sd_inv=parameters(end);
%stoch_inv_cost=sqrt(2)*repmat(sd_inv,1,N,1).*reshape(x_inv,1,1,n_node_inv);%1*N*n_node_inv

stoch_inv_cost=zeros(1,N,n_node_inv);

  [inv_cost]=inv_cost_func(k_t,I_t,stoch_inv_cost,theta);


%%% Value iteration
spec_V_iter=[];
spec_V_iter.ITER_MAX=OPI_param;
spec_V_iter.TOL=1e-9;
spec_V_iter.update_spec=spec.update_spec;
spec_V_iter.Anderson_acceleration=1;

if 1==0
    %%% Without precompute basis_t1
    [output,other_vars,iter_info]=spectral_func(...
        @V_update_func,spec_V_iter,...
        {V_t},...
        I_t,inv_cost,pi_mat,...
        k_t,exo_t,exo_t1,basis_t,inv_multiply_t,basis_exo_t1,...
        x_inv,w_inv,...
        state_min,state_max,Smol_elem,mu_max,d,ind);
else
     %%% Precompute basis_t1
     coef_approx_V=[];
     if spec_precompute==1
        [basis_t1]=...
            V_t1_func_precompute((1-delta_param)*reshape(k_t,n_pts,N,1)+I_t,...
            basis_exo_t1,coef_approx_V,state_min,state_max,Smol_elem,mu_max,ind,w_inv);

     else
        [basis_t1]=...
            V_t1_func((1-delta_param)*reshape(k_t,n_pts,N,1)+I_t,...
        exo_t1,coef_approx_V,state_min,state_max,Smol_elem,mu_max,ind,w_inv);
     end

     if 1==0
         [output,other_vars,iter_info]=spectral_func(...
            @V_update_func_given_basis_t1,spec_V_iter,...
            {V_t},...
            I_t,basis_t1,inv_cost,pi_mat,...
            k_t,exo_t,exo_t1,basis_t,inv_multiply_t,basis_exo_t1,...
            x_inv,w_inv,...
            state_min,state_max,Smol_elem,mu_max,d,ind);

     else
        % inv_cost:n_pts*N*n_node_inv
        E_inv_cost=reshape(sum(inv_cost.*reshape(w_inv,1,1,n_node_inv),3),n_pts,N);

        func_for_krylov_anonymous= @(V_vec)func_for_krylov(V_vec,...
            I_t,basis_t1,inv_cost,pi_mat,...
            k_t,exo_t,exo_t1,basis_t,inv_multiply_t,basis_exo_t1,...
            x_inv,w_inv,...
            state_min,state_max,Smol_elem,mu_max,d,ind);
        
        iter_info=[];
        ITER_MAX_gmres=100;
        TOL_gmres=1e-6;
        [V_new_vec,flag_vec,relres,iter_gmres,resvec] = gmres(func_for_krylov_anonymous,...
            pi_mat(:)-E_inv_cost(:),[],...
            TOL_gmres,ITER_MAX_gmres,[],[],V_t(:)); % solve for Krylov vector
        V_t_updated=reshape(V_new_vec,n_pts,N);
     end

end

if isempty(iter_info)==0
    V_t_updated=output{1};
    veval_total=veval_total+n_pts*N*iter_info.feval;
else
    veval_total=veval_total+n_pts*N*prod(iter_gmres);
end


coef_approx_V=inv_multiply_t*V_t_updated;
[I_t_updated,basis_t1,geval]=Newton_func0(I_t,k_t,exo_t1,basis_exo_t1,stoch_inv_cost,...
    coef_approx_V,state_min,state_max,Smol_elem,mu_max,d,ind,w_inv,theta,I_min,I_max);

geval_total=geval_total+geval;

out={I_t_updated,V_t_updated};
other_vars=[];

end
