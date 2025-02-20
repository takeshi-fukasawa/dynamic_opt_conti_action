function [] = ccprofit(nfirms,descn)
% competition
% local i, n, q, p;

  global ggamma D f
  global profit agprof csurplus share pmcmarg concent
  i = 1;
  while i <= descn;
    progress(i);
    w = cdecode(i,nfirms+1);
    theta = ggamma * exp(-w);  % marginal cost

%   quan = solveeq(theta); inlined
%   Solve for equilibrium with n firms; reduce n until all firms
%   want to produce quantity > 0

    n=nfirms;
    p = (D + sum(theta(1:n)))/(n+1);
    while ~((p - theta(n) >= 0) | (n==1));
      n=n-1;
      p = (D + sum(theta(1:n)))/(n+1);
      end
    q = zeros(nfirms,1);
    if p - theta(n) > 0;
      q(1:n) = p - theta(1:n);
      end
    quan=q;

    pstar = D - sum(quan);   % Equilibrium price
    profstar = (pstar>theta).*(pstar-theta).*quan - f; % Equilibrium profits

    profit(i,:) = profstar';
    csurplus(i) = 0.5*sum(quan)*sum(quan);
    agprof(i,:) = profstar';
    share(i,:) = quan';
    if sum(quan) > 0;
      pmcmarg(i) = pstar / (sum(theta.*quan)) * sum(quan);
      concent(i) = max(quan)/sum(quan);
    else;
      pmcmarg(i) = 1;
      end
    i = i+1;
    end
end
