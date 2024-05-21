% oct2mat
% gauss2oct

% This version is oneagentdescstat.prg
% Written by: Gautam Gowrisankaran
% April 21, 1993
% This program generates the output programs used for descriptive statistics,
% for the multiplant monopolist/social planner value functions


function [] = ds_sa()


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

disp(sprintf('\nENTRY-EXIT SIMULATION\n'));
disp(['  Periods to simulate:', sprintf(' %6d', numtimes), ...
  '        Initial state:', sprintf(' %2d', wstart) ]);
disp('  Initializing ...');

% In the social planner problem, k=0 is out; k=1 is the lowest state
% That is why kmax,entry_k is one higher than in the MPNE problem.

kmax = kmax+1;
entry_k = entry_k+1;
wstart = (wstart>0).*(wstart+1);

% Set up binomial coefficients for decoding/encoding of n-tuples

binom = eye(rlnfirms+kmax+1);
binom = [zeros(rlnfirms+kmax+1,1),binom];
i=2;
while i <= rlnfirms+kmax+1;
  binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
  i=i+1;
  end

wmax = binom(rlnfirms+kmax+1,kmax+2);

% Load in all the data stored by the equilibrium generation program
% This data is: x (investment), p (probability of state rising),
%   isentry, whichin



load(['a.' c.PREFIX '_oneag' int2str(rlnfirms) '.mat']);
x = newx; p = prising; isentry; whichin;

% Load in all data from the static profit calculation
% The data is: firm profits, consumer surplus, market shares at each state,
% price/cost margins, one-firm concentration ratios.



load(['a.' c.PREFIX '_cons' int2str(rlnfirms) '.mat'])
profit = agprof; csurplus; share; pmcmargm = pmcmarg; concentm = concent;

active = zeros(rlnfirms+1,1);  % No. of periods with n firms active.
exitors = 0; % No. of periods with exit
entrants = 0; % No. of periods with entry
entexit = 0;
invest = zeros(numtimes,1);  % Average investment
pmcmarg = zeros(numtimes,1); % Price/mc margin
concent = zeros(numtimes,1); % One-firm concentration ratio
lifedis = 0;  % Distribution of firm lifespans
valuedis = 0;  % Distribution of firm's total profits
lifemx = zeros(rlnfirms,1); % No. of periods that firm has been active for
valuemx = zeros(rlnfirms,1); % Total profits of each firm to date

wthis = wstart;
lastsize = share(encode(wthis),:)';  % Shares of firms last iteration
t = 0;
while t < numtimes;
  codew = encode(wthis);
  % Figure out exit

  i=whichin(codew);
  wtrans = zeros(rlnfirms,1);
  if i > 0;
    wtrans(1:i) = wthis(1:i);
    end

  y1 = sum(wtrans > zeros(rlnfirms,1));
  y2 = sum(wthis > zeros(rlnfirms,1));

  % Complete lifespan and value dists for firms that are exiting, if any.

  if y2 > y1;
    % "Exit in period " t;

    lifedis = [lifedis;lifemx(y1+1:y2)];
    valuemx(y1+1:y2) = valuemx(y1+1:y2)+(beta^(lifemx(y1+1:y2)+1))*phi;
    valuedis = [valuedis;valuemx(y1+1:y2)];

    lifemx(y1+1:y2) = zeros(y2-y1,1);
    valuemx(y1+1:y2) = zeros(y2-y1,1);
    end
  codew2 = encode(wtrans);
  prob = p(codew,:)';  % Probability of investment causing a rise in eff level

  lifemx = lifemx + (wtrans > 0);
  valuemx = valuemx + (beta.^lifemx).*(-x(codew2,:)' + profit(codew2,:)');

  % Now figure out entry

  entrypr = rand(1,1);
  yesentry = (isentry(codew2)>entrypr);
  if yesentry;
    % "Entry in period " t;

    % prob(rlnfirms) = 0;
    wtrans(rlnfirms) = entry_k;
    entryfee = x_entryl + entrypr * (x_entryh - x_entryl);
    valuemx(rlnfirms) = -entryfee;

  else; entryfee = 0;
    end

  wnext = wtrans+(prob>=rand(rlnfirms,1)) - (rand(1,1) <= delta);
  wnext = max(([wnext,zeros(rlnfirms,1)])')';

  % Now, tally the statistics

  active(y1+1) = active(y1+1)+1;
  exitors = exitors + (y2 > y1);
  entrants = entrants + yesentry;
  entexit = entexit + ((yesentry) & (y2 > y1));
  thissize = share(codew2,:)';
  l = sum(lastsize);
  concent(t+1) = concentm(codew2);
  invest(t+1) = sum(x(codew2,:)');
  pmcmarg(t+1) = pmcmargm(codew2);

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
% This procedure takes a weakly descending n-tuple (n = nfirms), with
% min. elt. 0, max. elt. kmax, and encodes it into an integer
% local code,digit,i;

  global rlnfirms binom
  code = 1; % Coding is from 1 to wmax

  i = 1;
  while i <= rlnfirms;
    digit = ntuple(i);
    code = code + binom(digit+rlnfirms-i+1,digit+1);
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
  while i<=nbreak;
    bp(i)=bp(i-1)+step;
    i=i+1;
    end
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