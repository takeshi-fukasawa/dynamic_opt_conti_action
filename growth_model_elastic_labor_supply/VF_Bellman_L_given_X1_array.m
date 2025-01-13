
function [out,other_vars] = VF_Bellman_L_given_X1_array(V,profit,X1_array,X0,beta,weight_nodes)

global relative_V_spec

[n_grid,n_coef,n_nodes]=size(X1_array);
vf_coef=X0\V;

%%% X1: n_grid*n_coef*n_nodes
%%% vf_coef: n_coef*1
EV=reshape(...
    sum(reshape(X1_array,n_grid,n_coef,n_nodes).*reshape(vf_coef,1,n_coef,1),2),...
    n_grid,n_nodes);

V_new = profit+beta*EV*weight_nodes; % Bellman equation        

if relative_V_spec>=1
    V_new=relative_V_func(V_new,relative_V_spec);
end

out={V_new};
other_vars=[];
end
