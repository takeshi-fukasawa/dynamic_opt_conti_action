function results=create_results_func(iter_info)

    DIST=iter_info.DIST_table(iter_info.feval,:);

    results=[round(iter_info.t_cpu,2),iter_info.feval,...
        iter_info.geval_total,...
        log10(DIST)];

end
