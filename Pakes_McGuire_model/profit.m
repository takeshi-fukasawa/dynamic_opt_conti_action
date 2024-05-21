% oct2mat
% gauss2oct


function [] = profit(c)


% function [] = cqprofit(nfirms, descn)
% competition
% local i;
% i = 1;
% while i <= descn;
%   progress(i);
%   w = qdecode(i,nfirms+1);
%   p = newton(p,&cfunk);
%   profstar = M*p.*sigma - M*mc*sigma;
%   profit(i,:) = profstar';
%   csurplus(i) = M*ln(1+sum(exp(ln(egw)-p)));
%   agprof(i,:) = profstar';
%   share(i,:) = sigma';
%   pmcmarg(i) = sum(p.*sigma) / mc / sum(sigma);
%   concent(i) = max(sigma)/sum(sigma);
%   i = i+1;
% end
% endfunction


% function [] = mqprofit(nfirms, descn)
% monopoly
% local i;
% i = 1;
% while i <= descn;
%   progress(i);
%   w = qdecode(i,nfirms+1);
%   w = max(([(w-3),(-7*ones(nfirms,1))])');
%     As zero here represents being out, move everything down by one,
%     except for zero, as there is no negative efficiency level
%   p = newton(p,&mfunk);
%   profstar = m*p.*sigma - m*mc*sigma;
%   profit(i) = sum(profstar);
%   agprof(i,:) = profstar';
%   csurplus(i) = M*ln(1+sum(exp(ln(egw)-p)));
%   share(i,:) = sigma';
%   pmcmarg(i) = sum(p.*sigma) / mc / sum(sigma);
%   concent(i) = max(sigma)/sum(sigma);
%   i = i+1;
% end
% endfunction


% function [] = sqprofit(nfirms, descn)
% social planner
%local i, denom, n;
% i = 1;
% while i <= descn;
%   progress(i);
%   w = qdecode(i,nfirms+1);
%   w = max(([(w-3),(-7*ones(nfirms,1))])');
%     As zero here represents being out, move everything down by one, except for
%     zero, as there is no negative efficiency level
%   p = mc*ones(nfirms,1);
%   egw = eg(w);
%   egwp = exp(ln(egw)-p);
%   n = egw.*exp(-p);
%   denom = 1 + sum(n);
%   sigma = n./denom;
%   profit(i) = M*ln(1+sum(egwp));  % Consumer surplus
%   share(i,:) = sigma';
%   pmcmarg(i) = sum(p.*sigma) / mc / sum(sigma);
%   concent(i) = max(sigma)/sum(sigma);
%   i = i+1;
% end
% csurplus=profit;
% endfunction


% function [] = cfunk(p)
%   used for quality competition profit function
%   local n,denom;
%   egw = eg(w);
%   n = egw.*exp(-p);
%   denom = 1.0 + sum(n);
%   sigma = n./denom;
%   retp(-(p-mc).*sigma.*(1-sigma) + sigma);
%   retp(-(p-mc).*(1-sigma) + 1);
%   endfunction


% function [out1] = mfunk(p)
% used for quality-investment, monopoly model

%   Calculate the profit derivative wrt price for the multi-plant monopolist.
%   The derivative is as in (18) in PGM paper. However, we can simplify this
%   expression as follows:
%   first, divide through by sigma(n).
%   left with: FOCn = -(pn-mc)(1-sigman) + 1 + sum(k<>n)(sigmak(pk-mc))
%   second, define A = sum(pk-mc).
%   left with: FOCn = -(pn-mc)(1-sigman) + 1 + A - sigman(pn-mc).
%   ==> FOCn = -(pn-mc) + 1 + A
%   ==> FOC = -(p-mc) + 1 + A.

%   local n,denom,A;
%   Might need to fix this
%   egw = eg(w);
%   n = egw.*exp(-p);
%   denom = 1.0 + sum(n);
%   sigma = n./denom;
%   A = sum(sigma.*(p-mc));
%   out1 = -(p-mc) + 1 + A;
%   retp((p-mc).*(-1 + sum(sigma)) + 1);
%   endfunction


% function [] = eg(w)
% used for quality competition profit function
% Calculates e^g(w)

%   local i,wret;
%   wret = zeros(rows(w),1);
%   i=1;
%   while i <= rows(w);
%     if w(i) <= wstar;
%       wret(i) = exp(w(i));
%     else; wret(i) = exp(wstar)*(2.0-exp(-(w(i)-wstar)));
%     end
%     i=i+1;
%   end
%   retp(wret);
%   endfunction








% function [] = cpprofit(nfirms,descn)
% competition
% local i, q, n, sub, cap;
%   i = 1;
%   while i <= descn;
%     progress(i);
%     w = pdecode(i,nfirms+1);
%     quan = solveeq(w*tau); inlined
%     cap=w*tau;
%     Solve for equilibrium with n firms; reduce n until all firms
%     want to produce quantity > 0

%     n=nfirms;
%     sub = 0;
%     q = D / (n + 1) * ones(n,1);
%     while q(n) > cap(n);  % Capacity constraint violated

%       q(n) = cap(n);
%       sub = sub + q(n);
%       Find reduced demand curve

%       n=n-1;
%       if n == 0; break; end
%       q(1:n) = (D - sub) / (n + 1) * ones(n,1);
%       end
%     quan=q;
%     pstar = D - sum(quan);   % Equilibrium price
%     profstar = pstar.*quan; % Equilibrium profits

%     profit(i,:) = profstar';
%     csurplus(i) = 0.5*sum(quan)*sum(quan);
%     agprof(i,:) = profstar';
%     share(i,:) = quan';
%     if sum(quan) > 0;
%       pmcmarg(i) = (pstar+mc)/mc;
%       concent(i) = max(quan)/sum(quan);
%     else;
%       pmcmarg(i) = 1;
%       end
%     i = i+1;
%     end
% endfunction


% function [] = mpprofit(nfirms, descn)
% monopoly
% local i, j, totprod, q, cap;
% Monopolist produces quantity D/2, and charges price on demand curve, unless
% total capacity is less than D/2, in which case quantity = capacity

% i = 1;
% while i <= descn;
%   progress(i);
%   w = pdecode(i,nfirms+1);
%   w = max(([(w-1),zeros(nfirms,1)])');
%   As zero here represents being out, move everything down by one, except for
%   zero, as there is no negative efficiency level

%   cap = w*tau;
%   quan = solveeq(cap); inlined

%   Solve for total production, and then allocate the production among the
%   firms putting the maximum possible production on the first firms.

%   q = zeros(nfirms,1);
%   totprod = min([sum(cap);0.5*D]);
%   j=1;
%   while (totprod > 0) and (j <= nfirms);
%     q(j) = min([totprod;cap(j)]);
%     totprod = totprod - q(j);
%     j=j+1;
%     end
%   quan=q;
%   pstar = D-sum(quan);   % One-plant monopolist price
%   profit(i) = sum(quan)*pstar;       % Monopolist profits

%   agprof(i,:) = pstar*quan';
%   pmcmarg(i) = (pstar + mc)/mc;
%   csurplus(i) = 0.5*sum(quan)*sum(quan);
%   if sum(quan) > 0;
%     concent(i) = quan(1)/sum(quan);
%     share(i,:) = quan';
%     end
%   i=i+1;
%   end
% endfunction


% function [] = spprofit(nfirms,descn)
% social planner
% local i, j, cap, q, totprod;
% Social planner produces quantity D, and charges price on demand curve, unless
% total capacity is less than D, in which case quantity = capacity

% i = 1;
% while i <= descn;
%   progress(i);
%   w = pdecode(i,nfirms+1);
%   w = max(([(w-1),zeros(nfirms,1)])');
%   As zero here represents being out, move everything down by one, except for
%   zero, as there is no negative efficiency level

%   cap = w*tau;
%   quan = solveeq(cap); inlined
%   Solve for total production, and then allocate the production among the
%   firms putting the maximum possible production on the first firms.

%   q = zeros(nfirms,1);
%   totprod = min([sum(cap);D]);
%   j=1;
%   while (totprod > 0) and (j <= nfirms);
%     q(j) = min([totprod;cap(j)]);
%     totprod = totprod - q(j);
%     j=j+1;
%     end
%   quan=q;
%   pstar = D-sum(quan);   % One-plant monopolist price

%   agprof(i,:) = pstar*quan';
%   csurplus(i) = 0.5*sum(quan)*sum(quan);
%   profit(i) = sum(agprof(i,:)')+csurplus(i); % Actually, total surplus

%   pmcmarg(i) = (pstar + mc)/mc;
%   if sum(quan) > 0;
%     concent(i) = quan(1)/sum(quan);
%     share(i,:) = quan';
%     end
%   i=i+1;
%   end
% endfunction


% function [out1] = puttuple(ntuple,place,val)
% This procedure puts val in position 'place' of ntuple, and then reorders
% ntuple to make sure that it is in descending order

%   ntuple(place) = val;
%   out1 = flipud(sortrows(flipud(ntuple),1));
%   endfunction


% function [out1] = qdecode(code,nfirms)
% Bertrand
% This procedure takes a previously encoded number, and decodes it into
% a weakly descending n-tuple (n = nfirms - 1)
% local ntuple,digit,i;

% code = code-1;
% ntuple = zeros(nfirms-1,1);
% i = 1;
% while i <= nfirms - 1;
%   digit = 0;
%   while binom(digit+nfirms-i+1,digit+2) <= code;
%     digit=digit+1;
%   end
%   ntuple(i) = digit;
%   code = code-binom(digit+nfirms-i,digit+1);
%   i = i+1;
% end

% Now convert to format of starting at -7, and jumping by 3's

% ntuple = (ntuple.*3 - 7);
% out1 = ntuple;
% endfunction






% function [out1] = newton(p,&objfunk)
% This procedure performs a simple Newton-Raphson search to find the root
% of the function objfunk
% local objfunk:proc,deriv,pnew,epsilon,norm,tol,x,i,id,iter,maxiter;

% id = eye(rows(p));
% epsilon=0.0001; tol = 1e-8;
% maxiter = 100; iter = 0;
% norm = tol+1;
% deriv = zeros(rows(p),rows(p));
% while (norm > tol) and (iter < maxiter);
%   iter=iter+1;
%   Calculate function at p

%   x = objfunk(p);
%   Calculate derivative matrix

%   i=1;
%   while i <= rows(p);
%     deriv(:,i) = (objfunk(p+(epsilon.*id(:,i)))-x)/epsilon;
%     i=i+1;
%   end
%   pnew = p - 0.6*inv(deriv)*x;
%   norm = max(abs(pnew - p));
%   p = pnew;
% end
% if norm > tol;
%   "Error: Newton method could not solve this problem.";
% end
% out1 = pnew;
% endfunction


% #include pmg.h
% This is a shell for different types of profit functions

global c binom p mc sigma wstar M
global ggamma D f
global profit agprof csurplus share pmcmarg concent

nfmax = c.MAX_FIRMS; % max # of active firms
kkmax = c.KMAX;      % max efficiency level attainable
it = c.IND_TYPE;     % investment type (quality/mc/capacity)
et = c.EQL_TYPE;     % equilibrium type (Nash/monopoly/social planner)

% "*** Computing profit function ***";
% Set up binomial coefficients for decoding/encoding of n-tuples

if strcmp(et, 'COMPETITION') kmax = kkmax;
else kmax = kkmax+1; end
% ? - might change this
if strcmp(it, 'CAPACITY') & strcmp(et, 'PLANNER') kmax = kmax+1; end

binom = eye(nfmax+kmax+2);
binom = [zeros(nfmax+kmax+2,1),binom];
i=2;
while i <= nfmax+kmax+2;
  binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
  i=i+1;
  end

nfirms=1;
while nfirms <= nfmax;
%     Number of descending n-tuples

  disp(sprintf('\nFirms: %d', nfirms))
  descn = binom(nfirms+kmax+1,kmax+2);
  disp(sprintf('Industry structures to compute: %d', descn))
  if strcmp(et, 'COMPETITION') nagents = nfirms;
  else nagents = 1; end

  profit = zeros(descn,nagents);
% The following variables are used for comparative statics
  agprof = zeros(descn,nfirms); % Aggregate profits per industry structure
  csurplus = zeros(descn,1);    % Consumer surplus
  share = zeros(descn,nfirms);  % Market share of each firm
  pmcmarg = zeros(descn,1);     % Price/mc margin, average by sales
  concent = zeros(descn,1);     % One-firm concentration ratio


%      now, call appropriate profit function

    if strcmp(it, 'QUALITY');
%      mc = c.MC;
%      M = c.MKT_SIZE;
%      wstar = c.WSTAR;
%      w = []; egw = []; egwp = []; p = []; profstar = [];
%      sigma=zeros(nfirms,1);
%      p=5.5*ones(nfirms,1);
%      if strcmp(et, 'COMPETITION');
%         cqprofit(nfirms, descn);
%      elseif strcmp(et, 'MONOPOLY');
%         mqprofit(nfirms,descn);
%      elseif strcmp(et, 'PLANNER');
%         sqprofit(nfirms,descn);
%      end
    elseif strcmp(it, 'COST');
       D = c.INTERCEPT;
       f = c.FIXED_COST;
       ggamma = c.GAMMA;
       quan = []; profstar = []; w = []; theta = []; pstar = [];
       if strcmp(et, 'COMPETITION');
         ccprofit(nfirms, descn);
       elseif strcmp(et, 'MONOPOLY');
         mcprofit(nfirms,descn);
       elseif strcmp(et, 'PLANNER');
         scprofit(nfirms,descn);
         end
    elseif strcmp(it, 'CAPACITY');
%      D = c.INTERCEPT;
%      mc = c.MC;
%      tau = c.TAU;
%      quan = []; profstar = []; w = []; pstar = [];
%      if strcmp(et, 'COMPETITION');
%         cpprofit(nfirms, descn);
%      elseif strcmp(et, 'MONOPOLY');
%         mpprofit(nfirms,descn);
%      elseif strcmp(et, 'PLANNER');
%         spprofit(nfirms,descn);
%      end
       end

%   write output
%     "Generating output --> ";;file1;;", ";;file2;

  s = int2str(nfirms);
  %save(['a.' c.PREFIX '_pr' s '.mat'], 'profit');
  %save(['a.' c.PREFIX '_cons' s '.mat'], ...
  %  'agprof', 'csurplus', 'share', 'pmcmarg', 'concent')

  nfirms=nfirms+1;
  end
c.PROFIT_DONE=1;



function [] = progress(i)
% report progress
  if mod(i, 50) == 0;
    % "Computed: ";;compact(i);;"\r";;
    disp(sprintf('  Computed: %d', i));
    end


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



function [out1] = cdecode(code,nfirms)
% Cournot
% This procedure takes a previously encoded number, and decodes it into
% a weakly descending n-tuple (n = nfirms - 1)
% local ntuple,digit,i;

  global binom
  code = code-1;
  ntuple = zeros(nfirms-1,1);
  i = 1;
  while i <= nfirms - 1;
    digit = 0;
    while binom(digit+nfirms-i+1,digit+2) <= code;
      digit=digit+1;
      end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i,digit+1);
    i = i+1;
    end

% Now convert to format of starting at -4, and jumping by 1's

  ntuple = ntuple-4;
  out1 = ntuple;


function [out1] = pdecode(code,nfirms)
% Capacity
% This procedure takes a previously encoded number, and decodes it into
% a weakly descending n-tuple (n = nfirms - 1)
% local ntuple,digit,i;

  global binom
  code = code-1;
  ntuple = zeros(nfirms-1,1);
  i = 1;
  while i <= nfirms - 1;
    digit = 0;
    while binom(digit+nfirms-i+1,digit+2) <= code;
      digit=digit+1;
    end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i,digit+1);
    i = i+1;
  end
  out1 = ntuple;
