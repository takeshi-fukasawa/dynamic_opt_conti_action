function [obj,obj_diff]=Q_func(x_j,v1,v2,a,beta,INV_COST,QUAD_INV_COST)

    [FOC,p]=FOC_func(x_j,v1,v2,a,beta,INV_COST,QUAD_INV_COST);

    Q_temp=- INV_COST.*x_j-QUAD_INV_COST.*x_j.^2 + beta*(v1*p + v2*(1-p));

    obj=-Q_temp;
    obj_diff=-FOC;

end