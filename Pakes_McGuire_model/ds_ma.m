% oct2mat
% gauss2oct

% Entry & exit statistics for multiple agents problem


% This version is markds.prg
% Written by: Gautam Gowrisankaran
% April 25, 1993
% This program generates the output programs used for the Markov-perfect game.


function [] = ds_ma()


% #include pmg.h;
% "**** Computing Entry and Exit Statistics ****";
% #include init.h;

global c rlnfirms binom
kmax = c.KMAX;
x_entryl = c.ENTRY_LOW;
x_entryh = c.ENTRY_HIGH;
phi = c.SCRAP_VAL;
entry_k = c.ENTRY_AT;
rlnfirms = c.MAX_FIRMS;
beta = c.BETA;
delta = c.DELTA;

wstart = zeros(rlnfirms,1);

% print "Enter initial efficiency level(s)";
% print "Valid range for parameter values is 0 - " compact(kmax);
% i=1;
% while i<=rlnfirms;
%   wdef=iif(i==1, entry_k+2, 0);
%   prompt="Firm "$+compact(i);
%   wstart(i)=getint(prompt, 0, kmax, wdef);
%   i=i+1;
%   end
% numtimes=getint("Number of periods to simulate (1-50K, default 10K)", ...
%   1,50000,10000);

wstart = c.DS_WSTART;
numtimes = c.DS_NSIMX;

% Set up binomial coefficients for decoding/encoding of n-tuples

binom = eye(rlnfirms+kmax+1);
binom = [zeros(rlnfirms+kmax+1,1),binom];
i=2;
while i <= rlnfirms+kmax+1;
  binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
  i=i+1;
  end

disp(sprintf('\nENTRY-EXIT SIMULATION\n'));
disp(['  Periods to simulate:', sprintf(' %6d', numtimes), ...
  '        Initial state:', sprintf(' %2d', wstart) ]);
disp('  Initializing ...');

wmax = binom(rlnfirms+1+kmax,kmax+2);

% Load in all the data stored by the equilibrium generation program
% This data is: v (value), x (investment), p (probability of state rising),
% isentry



load(['a.' c.PREFIX '_markov' int2str(rlnfirms) '.mat']);
v = newvalue; x = newx; p = prising; isentry;

% The data is: firm profits, consumer surplus, market shares at each state,
% price/cost margins, one-firm concentration ratios.



load(['a.' c.PREFIX '_cons' int2str(rlnfirms) '.mat'])
profit = agprof; csurplus; share; pmcmargm = pmcmarg; concentm = concent;

active = zeros(rlnfirms+1,1);  % No. of periods with n firms active.
exitors = 0; % No. of exitors
entrants = 0; % No. of entrants
entexit = 0; % No. of periods with both entry and exit
invest = zeros(numtimes,1);  % Average investment
pmcmarg = zeros(numtimes,1); % Price/mc margin
concent = zeros(numtimes,1); % One-firm concentration ratio, in each period
lifedis = 0;  % Distribution of firm lifespans
valuedis = 0;  % Distribution of firm's total profits
lifemx = zeros(rlnfirms,1); % No. of periods that firm has been active for
valuemx = zeros(rlnfirms,1); % Total profits of each firm to date

wthis = wstart;
lastsize = share(encode(wthis),:)';  % Shares of firms last iteration
t = 0;
while t < numtimes;
  % Get probabilities of investment causing a rise in eff level, as well
  % as actual investment and value function for this state tuple

  codew2 = encode(wthis);
  pr = p(codew2,:)';
  xx = x(codew2,:)';
  vv = v(codew2,:)';

  % Figure out exit

  wtrans = zeros(rlnfirms,1);
  i = (min(vv) == phi)*(minind(vv)-1) + (min(vv) > phi)*rlnfirms;
  if i > 0;
    wtrans(1:i) = wthis(1:i);
    end

  lifemx = lifemx + (wthis > 0);
  thisexit = (wtrans == 0) & (lifemx > 0);
  % if we want to force exit in the last period, replace by
  % thisexit = iif(t==numtimes-1,
  % (lifemx > 0),
  % (wtrans == 0) .and (lifemx > 0));

  lifemx = lifemx - (wthis > 0);

  % Complete lifespan and value dists for firms that are exiting, if any.

  if sum(thisexit) > 0;
    %"Exit in period " t;

    % lifedis = [lifedis;selif(lifemx,thisexit)];
    % valval = selif(valuemx,thisexit)+(beta^(1+selif(lifemx,thisexit)))*phi;
    [ii,jj,v1] = find(lifemx.*thisexit);
    [ii,jj,v2] = find(valuemx.*thisexit);
    lifedis = [lifedis; v1];
    valval = v2+(beta.^(1+v1))*phi;

    valuedis = [valuedis;valval];
    lifemx = lifemx .* (1-thisexit);
    valuemx = valuemx .* (1-thisexit);
    end

  codew = encode(wtrans);
  valuemx = valuemx + (beta.^lifemx).*(-xx + profit(codew,:)');
  numin = sum(wtrans>0);

  % Now figure out entry

  yesentry = 0;
  entrypr = rand(1,1);
  yesentry = (isentry(codew)>entrypr);
  if yesentry;
    % pr(rlnfirms) = 0;
    wtrans(rlnfirms) = entry_k;
    entryfee = x_entryl + entrypr * (x_entryh - x_entryl);
    valuemx(rlnfirms) = -entryfee;
    % "Entry in period " t "with fee" entryfee;
  else; entryfee = 0;
    end

  lifemx = lifemx + (wtrans > 0);
  wnext = wtrans+(pr>=rand(rlnfirms,1)) - (rand(1,1) <= delta);
  wnext = max(([wnext,zeros(rlnfirms,1)])')';

  % Now, tally the statistics

  active(numin+1) = active(numin+1)+1;
  exitors = exitors + sum(thisexit);
  entrants = entrants + yesentry;
  entexit = entexit + ((yesentry) & (sum(thisexit)));
  thissize = share(codew,:)';

  % Avoid division by zero for job creation %

  l = sum(lastsize);
  invest(t+1) = sum(xx);
  pmcmarg(t+1) = pmcmargm(codew);
  concent(t+1) = concentm(codew);

  if mod(t+1, 1000) == 0;
    disp([sprintf('  Periods simulated:   %6d        Current state:', t+1) ...
      sprintf(' %2d', wthis) ]);
    end

  % Now re-sort all firm level data, to reflect the fact that firms
  % must be in descending order next year

  temp = flipud(sortrows([wnext,lifemx,valuemx,thissize],1));
  wthis = temp(:,1); lifemx = temp(:,2);
  valuemx = temp(:,3); lastsize = temp(:,4);
  t = t+1;
  end

ii = [0:rlnfirms];
disp(sprintf('\nIndustry characterization\n'));
disp(sprintf('  Periods w/ %d firms active: %6d\n', [ii; active']));
disp(sprintf('  Periods w/ exit:           %6d', exitors));
disp(sprintf('  Periods w/ entry:          %6d', entrants));
disp(sprintf('  Periods w/ entry & exit:   %6d', entexit));

disp(sprintf('\n  Mean investment:     %8.2f (%8.2f)', ...
  mean(invest), std(invest)));
disp(sprintf('  Mean p-c margin:     %8.2f (%8.2f)', ...
  mean(pmcmarg), std(pmcmarg)));
disp(sprintf('  Mean 1-firm concent: %8.2f (%8.2f)', ...
  mean(concent), std(concent)));


if rows(lifedis) > 1;
  lifedis = lifedis(2:rows(lifedis));
  valuedis = valuedis(2:rows(valuedis));
  % (flipud(sortrows((lifedis~valuedis),2)));

  disp(sprintf('  Mean value:          %8.2f (%8.2f)', ...
    mean(valuedis), std(valuedis)));
  disp(sprintf('  Mean lifespan:       %8.2f (%8.2f)\n', ...
    mean(lifedis), std(lifedis)));
  disp(sprintf('  Total firms in history:    %6d', rows(lifedis)));

  end

temp = flipud(sortrows([lifemx,valuemx],1));
lifemx = temp(:,1);
valuemx = temp(:,2);
i = (min(lifemx) == 0)*(minind(lifemx)-1) + (min(lifemx) > 0)*rlnfirms;
if i > 0;
  disp([ sprintf('  Currently active firms have lived and earned:\n') ...
    sprintf('  %6d %8.2f\n', [lifemx(1:i)'; valuemx(1:i)']) ]);
  end

disp('  Lifespan Distribution');
distr(lifedis,20);

disp(sprintf('\n  Value Distribution'));
distr(valuedis,20);



function [out1] = encode(ntuple)
% This procedure takes a weakly descending n-tuple (n = rlnfirms), with
% min. elt. 0, max. elt. kmax, and encodes it into an integer
% local code,digit,i;

  global rlnfirms binom
  code = 1; % Coding is from 1 to wmax

  i = 1;
  while i <= rlnfirms;
    digit = ntuple(i);
    code = code + binom(digit+rlnfirms+1-i,digit+1);
    i=i+1;
  end
  out1 = code;


function [] = distr(x, nbreak)
  % local ncount, xmax, xmin, step, bp, i, j, total;

  xmax=max(x);
  xmin=min(x);
  bp=zeros(nbreak,1);
  step=(xmax-xmin)/nbreak;
  ncount=zeros(nbreak,1);
  bp(1)=xmin+step;
  i=2;
  while i<nbreak;
    bp(i)=bp(i-1)+step;
    i=i+1;
    end
  bp(nbreak)=xmax;
  i=1;
  while i<=rows(x);
    j=1;
    while bp(j)<x(i);
      j=j+1;
    end
    ncount(j)=ncount(j)+1;
    i=i+1;
    end
  total=sum(ncount);

  disp(sprintf('  %6.0f - %6.0f: %5.3f', xmin, bp(1), ncount(1)/total));
  j=2;
  while j<=nbreak;
    disp(sprintf('  %6.0f - %6.0f: %5.3f', bp(j-1), bp(j), ncount(j)/total));
    j=j+1;
    end