function [] = mcprofit(nfirms,descn)
% monopolist
% local i, numin;

  global ggamma D f
  global profit agprof csurplus share pmcmarg concent
  i = 1;
  while i <= descn;
    progress(i);

%   The monopolist will always choose to produce everything from the lowest
%   priced firm, so it acts like a 1-plant firm for the static profits.

    w = cdecode(i,nfirms+1);
    numin = sum(w > -4); % No. of firms in, for fixed-fee computation
    agprof(i,:) = -f * ((w > -4)');
    w = max(([(w-1),(-4*ones(nfirms,1))])');

%   As zero here represents being out, move everything down by one, except for
%   zero, as there is no negative efficiency level

    theta = ggamma * exp(-w(1));  % marginal cost
    pstar = 0.5*(D + theta);   % One-plant monopolist price
    quan = (pstar>theta)*(pstar-theta);
    profstar = quan*(pstar-theta) - f*numin; % Monopolist profits

    profit(i) = profstar;
    agprof(i,1) = agprof(i,1) + profstar;
    csurplus(i) = 0.5*quan*quan;
    pmcmarg(i) = pstar / theta;
    concent(i) = 1;
    if nfirms > 1;
      quan = [quan;zeros(nfirms-1,1)];
      end
    share(i,:) = quan';
    i = i+1;
    end

end
