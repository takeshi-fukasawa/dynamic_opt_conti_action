function results=create_results_func(iter_info)

    global profit
    n_grid=size(profit,1);

    DIST=iter_info.DIST_table(iter_info.feval,:);

    results=[round(iter_info.t_cpu,2),iter_info.feval,...
        iter_info.feval*n_grid,...
        iter_info.geval_total,...
        log10(DIST)];

end
