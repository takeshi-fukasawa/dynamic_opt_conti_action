
function [out,other_vars] = VF_Bellman_L_given_X1_array(V,profit,X1_array,X0,beta,weight_nodes)

[n_grid,n_coef,n_nodes]=size(X1_array);
vf_coef=X0\V;

%%% X1: n_grid*n_coef*n_nodes
%%% vf_coef: n_coef*1
EV=reshape(...
    sum(reshape(X1_array,n_grid,n_coef,n_nodes).*reshape(vf_coef,1,n_coef,1),2),...
    n_grid,n_nodes);

V_new = profit+beta*EV*weight_nodes; % Bellman equation        

out={V_new};
other_vars=[];
end
