
function [out,other_vars] = VF_Bellman_L_given_X1(V,profit,X1,X0,beta,weight_nodes)

global relative_V_spec

n_grid=size(V,1);
n_nodes=size(weight_nodes,1);
vf_coef=X0\V;

%%% X1: (n_grid*n_nodes)*n_coef
%%% vf_coef: n_coef*1
EV=X1*vf_coef;%(n_grid*n_nodes)*1

V_new = profit+beta*reshape(EV,n_grid,n_nodes)*weight_nodes; % Bellman equation        

if relative_V_spec>=1
    V_new=relative_V_func(V_new,relative_V_spec);
end

out={V_new};
other_vars=[];
end
