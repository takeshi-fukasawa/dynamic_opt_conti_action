function [out] =func_for_krylov(V,X1_array,X0,beta,weight_nodes)

    [n_grid,n_coef,n_nodes]=size(X1_array);
    vf_coef=X0\V;
    
    %%% X1: n_grid*n_coef*n_nodes
    %%% vf_coef: n_coef*1
    EV=reshape(...
        sum(reshape(X1_array,n_grid,n_coef,n_nodes).*...
        reshape(vf_coef,1,n_coef,1),2),...
        n_grid,n_nodes);
    
    out=V-beta*EV*weight_nodes;    
end
