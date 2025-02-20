% oct2mat
% gauss2oct

% This version is markwds.prg
% Written by: Gautam Gowrisankaran
% April 25, 1993
% This program generates the welfare output programs used for the
% Markov-perfect game.


function [] = welf_ma()


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

% Set up binomial coefficients for decoding/encoding of n-tuples

binom = eye(rlnfirms+kmax+1);
binom = [zeros(rlnfirms+kmax+1,1),binom];
i=2;
while i <= rlnfirms+kmax+1;
  binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
  i=i+1;
  end

% Number of possible industry structures

wmax = binom(rlnfirms+kmax+1,kmax+2);

% Load in all the data stored by the equilibrium generation program
% This data is: v (value), x (investment), p (probability of state rising),
%   isentry



load(['a.' c.PREFIX '_markov' int2str(rlnfirms) '.mat']);
v = newvalue; x = newx; p = prising; isentry;

% Load in all data from the static profit calculation
% The data is: firm profits, consumer surplus, market shares at each state,
% price/cost margins, one-firm concentration ratios.



load(['a.' c.PREFIX '_cons' int2str(rlnfirms) '.mat'])
profit = agprof; csurplus; share; pmcmargm = pmcmarg; concentm = concent;

consurp = zeros(numruns,1); % (Total) consumer surplus
prodsurp = zeros(numruns,1); % Producer surplus
totsurp = zeros(numruns,1);  % Total surplus

nr = 1;
while nr <= numruns;
  lifemx = zeros(rlnfirms,1);
  wthis = zeros(rlnfirms,1);
  wthis = wstart;
  t = 0;
  while t < numtimes;
    % Get probabilities of investment causing a rise in eff level, as well
    % as actual investment and value function for this state tuple

    codew2 = encode(wthis);
    pp = p(codew2,:)';
    xx = x(codew2,:)';
    vv = v(codew2,:)';

    % Find out which firms want to leave.

    wtrans = zeros(rlnfirms,1);
    i = (min(vv) == phi)*(minind(vv)-1) + (min(vv) > phi)*rlnfirms;
    if i > 0;
      wtrans(1:i) = wthis(1:i);
      end

    % Now figure out exit. Must capture people who exit voluntarily,
    % as well as firms whose efficiency has been driven down to zero.

    lifemx = lifemx + (wthis > 0);
    thisexit = (wtrans == 0) & (lifemx > 0);
    lifemx = lifemx - (wthis > 0);
    if sum(thisexit) > 0;
      lifemx = lifemx .* (1-thisexit);
      end
    codew = encode(wtrans);

    % Now figure out entry

    yesentry = 0;
    i=encode(wtrans);
    entrypr = rand(1,1);
    yesentry = (isentry(i)>entrypr);
    if yesentry;
      % pp(rlnfirms) = 0;
      wtrans(rlnfirms) = entry_k;
      entryfee = x_entryl + entrypr * (x_entryh - x_entryl);
    else; entryfee = 0;
      end
    wnext = wtrans+(pp>=rand(rlnfirms,1)) - (rand(1,1) <= delta);
    wnext = max(([wnext,zeros(rlnfirms,1)])')';

    lifemx = lifemx + (wtrans > 0);

    % Now, tally the statistics

    consurp(nr) = consurp(nr) + (beta^t)*csurplus(codew);
    prodsurp(nr) = prodsurp(nr) + ...
      (beta^t)*(sum(profit(codew,:)') ... % Static profits
      - entryfee + phi*sum(thisexit) ...  % Entry and exit fees
      - sum(xx))         ;                % Investment cost

    % Now re-sort all firm level data, to reflect the fact that firms
    % must be in descending order next year

    temp = flipud(sortrows([wnext,lifemx],1));
    wthis = temp(:,1); lifemx = temp(:,2);
    t = t+1;
    end

  if mod(nr, 10) == 0;
    disp([sprintf('  Runs simulated:      %6d', nr) ]);
    end
  nr=nr+1;
  end

totsurp = prodsurp + consurp;

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
