function [] = scprofit(nfirms, descn)
% local i, numin;

  global ggamma D f
  global profit agprof csurplus share pmcmarg concent
  i = 1;
  while i <= descn;
    progress(i);

%   The social planner will always choose to produce everything from the lowest
%   priced firm, so it acts like a 1-plant firm for the static profits.

    w = cdecode(i,nfirms+1);
    agprof(i,:) = -f*((w>-4)');
    numin = sum(w > -4); % No. of firms in, for fixed-fee computation
    w = max(([(w-1),(-4*ones(nfirms,1))])');

%   As zero here represents being out, move everything down by one, except for
%   zero, as there is no negative efficiency level

    theta = ggamma * exp(-w(1));  % marginal cost
    pstar = theta;  % Set price = mc, for social planner solution
    quan = (D>theta)*(D-pstar);
    profstar = 0.5*(D-theta)*quan; % Consumer surplus

    profit(i) = profstar - f*numin;  % Producer surplus
    csurplus(i) = profstar;
    pmcmarg(i) = 1;
    concent(i) = 1;
    if nfirms > 1;
      quan = [quan;zeros(nfirms-1,1)];
      end
    share(i,:) = quan';
    i = i+1;
    end
end
