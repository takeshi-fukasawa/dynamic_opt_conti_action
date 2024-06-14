% oct2mat
% gauss2oct


function [] = profit_func(c)

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
      mc = c.MC;
      M = c.MKT_SIZE;
      wstar = c.WSTAR;
      w = []; egw = []; egwp = []; p = []; profstar = [];
      sigma=zeros(nfirms,1);
      p=5.5*ones(nfirms,1);
      if strcmp(et, 'COMPETITION');
         cqprofit(nfirms, descn);
      elseif strcmp(et, 'MONOPOLY');
         mqprofit(nfirms,descn);
      elseif strcmp(et, 'PLANNER');
         sqprofit(nfirms,descn);
      end
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
  save(['./middle_output/' 'a.' c.PREFIX '_pr' s '.mat'], 'profit');
  save(['./middle_output/' 'a.' c.PREFIX '_cons' s '.mat'], ...
    'agprof', 'csurplus', 'share', 'pmcmarg', 'concent')

  nfirms=nfirms+1;
  end
c.PROFIT_DONE=1;

end
