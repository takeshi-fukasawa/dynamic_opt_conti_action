function [out] =func_for_krylov(V,X1,X0,beta,weight_nodes)

   global relative_V_spec
 
    n_grid=size(V,1);
    n_nodes=size(weight_nodes,1);
    vf_coef=X0\V;
    %%%vf_coef=(X0'*X0+0*eye(size(X0,2)))\(X0'*V);
    

    %%% X1: (n_grid*n_nodes)*n_coef
    %%% vf_coef: n_coef*1
    EV=X1*vf_coef;%(n_grid*n_nodes)*1
    EV=reshape(EV,n_grid,n_nodes);    

    if relative_V_spec==0
       out=V-beta*EV*weight_nodes;    
     else%relative_V_spec>=1
       out=V-relative_V_func(beta*EV*weight_nodes,relative_V_spec);
     end
end
