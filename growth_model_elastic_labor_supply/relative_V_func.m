function [vals_out]=relative_V_func(vals,relative_V_spec)

    global n_gridk n_grida

    if relative_V_spec==1
        vals_out=vals-vals(1);
    elseif relative_V_spec==2
        vals_out=reshape(vals,n_gridk,n_grida);
        vals_out=vals_out-mean(vals_out,1);
        vals_out=vals_out(:);
    end

end
