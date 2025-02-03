function [out] =func_for_krylov(V,X1_array,X0,beta,weight_nodes)

   global relative_V_spec
 
    [n_grid,n_coef,n_nodes]=size(X1_array);
    vf_coef=X0\V;
    
    %%% X1: n_grid*n_coef*n_nodes
    %%% vf_coef: n_coef*1
    EV=reshape(...
        sum(reshape(X1_array,n_grid,n_coef,n_nodes).*...
        reshape(vf_coef,1,n_coef,1),2),...
        n_grid,n_nodes);
    

    if relative_V_spec==0
       out=V-beta*EV*weight_nodes;    
     else%relative_V_spec>=1
       out=V-relative_V_func(beta*EV*weight_nodes,relative_V_func);
     end
end
