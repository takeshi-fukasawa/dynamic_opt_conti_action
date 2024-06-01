function out=iter_results_output_func(iter_info,resid_mat)
    resid_I=abs(reshape(resid_mat(:,1:2),[],1));
    resid_V=abs(reshape(resid_mat(:,3:4),[],1));
    
    out=[iter_info.t_cpu,...
    log10(mean(resid_I)),log10(max(resid_I)),log10(mean(resid_V)),log10(max(resid_V)),...
    iter_info.feval,1-iter_info.FLAG_ERROR];

end
