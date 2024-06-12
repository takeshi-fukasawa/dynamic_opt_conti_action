function [FOC,p]=FOC_func(x_j,v1,v2,a,beta,INV_COST,QUAD_INV_COST)

    p=(a.*x_j)./(1+a.*x_j);
    FOC=beta*a*((1-p).^2).*(v1-v2)-INV_COST-QUAD_INV_COST.*x_j;

end
