% oct2mat
% gauss2oct

% This version is oneagentwelfdescstat.prg
% Written by: Gautam Gowrisankaran
% April 21, 1993
% This program generates the output programs used for welfare descriptive


function [] = welf_sa()


% statistics, for the multiplant monopolist/social planner value functions

% #include pmg.h;
% "**** Computing Welfare Statistics ****";
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
% numtimes=getint("Number of periods to simulate (1-1K, default 100)", ...
%   1,1000,100);
% numruns=getint("Number of runs to simulate for (1-1K, default 100)", ...
%   1,1000,100);

wstart = c.DS_WSTART;
numtimes = c.DS_NSIMW;
numruns= c.DS_NRUNW;

jobstart = 5;  % Start of job creation statistics

disp(sprintf('\nWELFARE SIMULATION\n'));
disp(['  Number of periods:  ', sprintf(' %6d', numtimes), ...
  '        Initial state:', sprintf(' %2d', wstart) ]);
disp(['  Number of runs:     ', sprintf(' %6d', numruns) ]);
disp('  Initializing ...');

% In the social planner problem, k=0 is out; k=1 is the lowest state
% That is why kmax,entry_k is one higher than in the MPNE problem.

kmax = kmax+1;
entry_k = entry_k+1;
wstart = wstart+1;

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
%   price/cost margins, one-firm concentration ratios.



load(['a.' c.PREFIX '_cons' int2str(rlnfirms) '.mat'])
profit = agprof; csurplus; share; pmcmargm = pmcmarg; concentm = concent;

consurp = zeros(numruns,1); % (Total) consumer surplus
prodsurp = zeros(numruns,1); % Producer surplus
totsurp = zeros(numruns,1);  % Total surplus

nr = 1;
while nr <= numruns;
  wthis = zeros(rlnfirms,1);
  wthis = wstart;
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

    codew2 = encode(wtrans);

    % Probability of investment causing a rise in eff level

    prob = p(codew2,:)';

    % Now figure out entry

    entrypr = rand(1,1);
    yesentry = (isentry(codew2)>entrypr);
    if yesentry;
      % prob(rlnfirms) = 0;
      wtrans(rlnfirms) = entry_k;
      entryfee = x_entryl + entrypr * (x_entryh - x_entryl);
    else; entryfee = 0;
      end

    wnext = wtrans+(prob>=rand(rlnfirms,1)) - (rand(1,1) <= delta);
    wnext = max(([wnext,zeros(rlnfirms,1)])')';
    wnext = min(([wnext,kmax*ones(rlnfirms,1)])')';

    % Now, tally the statistics

    consurp(nr) = consurp(nr) + (beta^t)*csurplus(codew2);
    prodsurp(nr) = prodsurp(nr) + ...
      (beta^t)*(sum(profit(codew2,:)') ...    % Static profits
      - entryfee + phi*(y2 > y1)*(y2 - y1) ...% Entry and exit fees
      - sum(x(codew2,:)'));                   % Investment cost

    % Now re-sort all firm level data, to reflect the fact that firms
    % must be in descending order next year

    temp = flipud(sortrows(wnext,1));
    wthis = temp(:,1);
    t = t+1;
    end

  if mod(nr, 10) == 0;
    disp([sprintf('  Runs simulated:      %6d', nr) ]);
    end
  nr=nr+1;
  end

totsurp = consurp + prodsurp;

disp(sprintf('\nIndustry characterization\n'));
disp(sprintf('  Mean cons surplus:   %8.2f (%8.2f)', ...
  mean(consurp), std(consurp)));
disp(sprintf('  Mean prod surplus:   %8.2f (%8.2f)', ...
  mean(prodsurp), std(prodsurp)));
disp(sprintf('  Mean total surplus:  %8.2f (%8.2f)', ...
  mean(totsurp), std(totsurp)));

outsur = flipud(sortrows([consurp,prodsurp,totsurp],3));
disp([ sprintf('\nSurplus Realizations\n        CS       PS       TS\n') ...
  sprintf('  %8.2f %8.2f %8.2f\n', outsur') ]);



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
