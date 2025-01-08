function V_out=relative_V_func(V,relative_V_spec)

    global n_gridk n_grida

    if relative_V_spec==1
        V_out=V-V(1);
    elseif relative_V_spec==2
        V_out=reshape(V,n_gridk,n_grida);
        V_out=V_out-mean(V_out,1);
        V_out=V_out(:);
    end

end
