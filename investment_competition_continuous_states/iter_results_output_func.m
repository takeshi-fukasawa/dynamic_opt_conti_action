function out=iter_results_output_func(iter_info,resid_mat)
N=size(resid_mat,2)/2;
    resid_I=abs(reshape(resid_mat(:,1:N),[],1));
    resid_V=abs(reshape(resid_mat(:,N+1:2*N),[],1));
    
    out=[iter_info.t_cpu,...
    log10(mean(resid_I)),log10(max(resid_I)),...
    iter_info.feval,...
    iter_info.t_cpu./iter_info.feval,...
    iter_info.veval_total,...
    iter_info.geval_total,...
    1-iter_info.FLAG_ERROR];

end
