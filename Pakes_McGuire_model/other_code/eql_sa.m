% oct2mat
% gauss2oct

% This program is oneagent.prg. This calculates the optimal social planner
% or monopolist solution to a multi-plant industry.
% Written by: Gautam Gowrisankaran
% April 21, 1993



function [] = eql_sa()















% #include pmg.h;
% "**** Computing Dynamic Equilibrium ****";
% #include init.h;
% constants not modifiable by user

global c
global rlnfirms binom wmax encfirm multfac1 multfac2 etable1 etable2 oldx mask
global two_n kmax a delta oldvalue beta entry_k oneton x_entryl x_entryh profit
global phi isentry whichin
kmax = c.KMAX;
x_entryl = c.ENTRY_LOW;
x_entryh = c.ENTRY_HIGH;
phi = c.SCRAP_VAL;
entry_k = c.ENTRY_AT;
rlnfirms = c.MAX_FIRMS;
beta = c.BETA;
delta = c.DELTA;
a = c.INV_MULT;

tol = 0.1;  % Tolerance for convergence
kmax=kmax+1;
entry_k=entry_k+1;
newvalue = []; newx = [];

% For entry; no investment allowed for new plant

% oneton = seqa(1,1,rlnfirms);
oneton = [1:rlnfirms]';

% Set up binomial coefficients for decoding/encoding of n-tuples

binom = eye(rlnfirms+kmax+1);
binom = [zeros(rlnfirms+kmax+1,1),binom];
i=2;
while i <= rlnfirms+kmax+1;
  binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
  i=i+1;
  end

wmax = binom(rlnfirms+kmax+1,kmax+2);
disp(sprintf('\nFirms: %d   States: %d\nInitialization ...', rlnfirms, wmax));

% Isentry and whichin detail the probability of entry at any given state,
% and how plants the planner will have remain in production. These are not used
% by this program, but are used in generating comparative statics.

isentry = zeros(wmax,1); whichin = zeros(wmax,1);
encfirm = 3;  % Max. number of firms to encode in table
oldvalue = phi*ones(wmax,1); % The value function

i=1;
while i <= wmax;
  oldvalue(i,:) = oldvalue(i,:) + 0.1*i;
  i=i+1;
  end
oldx = zeros(wmax,rlnfirms); % Investment


load(['a.' c.PREFIX '_pr' int2str(rlnfirms) '.mat'])

two_n = 2^(rlnfirms-1);
dtable = [];
if rlnfirms > 1;
  % Build a mask of all binary numbers from 0 to two_n - 1

  mask = zeros(rlnfirms-1,two_n);
  msk = zeros(rlnfirms-1,1);
  i=1;
  while i <= two_n;
    mask(:,i)=msk;
    msk(rlnfirms-1)=msk(rlnfirms-1)+1;
    j=rlnfirms-1;
    while (msk(j) > 1) & (j > 1);
      msk(j) = 0; j=j-1; msk(j)=msk(j)+1;
      end
    i=i+1;
    end
  % print "Mask is " mask;
  end

% Make a table for quick decoding

dtable = zeros(rlnfirms,wmax);
i=1;
while i <= wmax;
  dtable(:,i) = decode(i);
   i=i+1;
  end

% Make tables for quick encoding

tempfirm = rlnfirms;
rlnfirms = min([rlnfirms;encfirm]);
multfac1 = (kmax+1).^(oneton(1:rlnfirms)-1);
if tempfirm > encfirm;
  multfac1 = [zeros(tempfirm-encfirm,1);multfac1];
  end
% Encode all numbers from 1 to kmax^rlnfirms

etable1 = zeros((kmax+1)^rlnfirms,1);
msk = zeros(rlnfirms,1);
i=1;
while i <= rows(etable1);
  etable1(i) = encode(flipud(sortrows(msk,1)));
  msk(rlnfirms)=msk(rlnfirms)+1;
  j=rlnfirms;
  while (msk(j) > kmax) & (j > 1);
    msk(j) = 0; j=j-1; msk(j)=msk(j)+1;
    end
  i=i+1;
  end

rlnfirms = tempfirm;

if rlnfirms > encfirm;
  multfac2 = [(kmax+1)^(oneton(1:rlnfirms-encfirm)-1);zeros(encfirm,1)];
  etable2 = zeros((kmax+1)^(rlnfirms-encfirm),1);
  msk=zeros(rlnfirms-encfirm,1);
  i=1;
  while i <= rows(etable2);
    etable2(i) = encode(flipud(sortrows(msk,1)));
    msk(rlnfirms)=msk(rlnfirms)+1;
    j=rlnfirms;
    while (msk(j) > kmax) & (j > 1);
      msk(j) = 0; j=j-1; msk(j)=msk(j)+1;
      end
    i=i+1;
    end
  etable2(i+1) = encode(flipud(sortrows([msk(1:rlnfirms-encfirm); ...
    zeros(encfirm,1)],1)))-1;
  end

disp(sprintf('Contraction ...'));
ix = 1;

norm = tol + 1.0;
while norm > tol;
  [newx,newvalue] = contract;
  norm = max(abs(oldvalue - newvalue));
  avgnorm = mean(abs(oldvalue-newvalue));
  oldx = newx; oldvalue = newvalue;
  disp(sprintf('  %2d    Sup norm: %8.4f      Mean norm: %8.4f', ...
    ix, norm, avgnorm));
  ix = ix+1;
  end

% Store data in file for inspection
% "Element/Newvalue (w): "; dtable'~newvalue; "";
% "Element/Newx (w x nfirms): "; dtable'~newx; "";
% "Element/Prob of p rising (w x nfirms): "; dtable'~(a.*newx./(1+a.*newx)); "";
% "Element/Isentry (w): "; dtable'~isentry; "";
% "Element/Whichin (w): "; dtable'~whichin; "";
% Store data in file for comparative statics program to read

prising = a.*newx./(1+a.*newx);
save(['a.' c.PREFIX '_oneag' int2str(rlnfirms) '.mat'], ...
   'newx', 'prising', 'isentry', 'whichin')
c.EQL_DONE = 1;



function [out1,out2] = contract()
  % This procedure does one iterative step on investment and the value fn
  % Implicit parameters are oldx, oldvalue
  % local v1,x1,i,w,locw,ww,www,
  % numin,   % Number of firms active for this w
  % locwmax,  % Wmax for the no. of firms being considered
  % tempfirm,
  % newx,newvalue,
  % probent;

  global wmax rlnfirms kmax binom phi isentry whichin
  newx = zeros(wmax,rlnfirms);
  newvalue = zeros(wmax,1);
  tempfirm = rlnfirms;
  locw = zeros(rlnfirms,1);

  % First figure out the value of all firms being inactive.

  ww=qencode(locw);
  [v1,x1,probent] = optimize(locw,ww);

  newvalue(ww) = v1;
  newx(ww,:) = x1';
  isentry(ww) = probent;
  whichin(ww) = 0;

  numin=1;
  while numin <= rlnfirms;
    locwmax = binom(numin+kmax+1,kmax+2);
    w=1;
    locw = zeros(rlnfirms,1);
    while w <= locwmax;
      rlnfirms = numin;
      locw(1:numin) = decode(w);
      rlnfirms = tempfirm;
      if locw(numin) == 0;
        w=w+1;
      else

        ww=qencode(locw);
        [v1,x1,probent] = optimize(locw,ww);

        % Now compare best value of having all the firms in, with the best
        % value of having all but one firm in.
        % On the firm we are at, all numin firms are in;
        % thus we can add phi to the
        % value of having one be out.

        locw(numin) = 0;
        isentry(ww) = probent;
        www = qencode(locw);
        newx(ww,:) = x1';
        if newvalue(www) + phi > v1;
          % Take the last firm out - this is best

          newvalue(ww) = newvalue(www)+phi;
          whichin(ww) = whichin(www);
        else;
          newvalue(ww) = v1;
          whichin(ww) = numin;
          end
        end

      w=w+1;
      end


    numin=numin+1;
    end

  out1 = newx;
  out2 = newvalue;


function [out1,out2,out3] = optimize(locw,ww)
% Zerolist contains new plants, so no investment allowed
% This procedure calculates optimal investment, and value fn., for a
% given  w (omega) combination
% Implicit parameters are oldx, oldvalue
% local newx, newvalue,i,firm,r,p,v1,v2,xothers,srtordr,locentr,
%   tempv1,tempv2,locwe,newxe,
%   newplace;  % Place of all w's after reordering for entry

  global oldx rlnfirms beta entry_k oneton x_entryl x_entryh profit a
  xothers = oldx(ww,:)';

  % First decide if entry is profitable

  if locw(rlnfirms) == 0;
    % Entry is possible

    [tempv1,v1] = calcval(locw,xothers,rlnfirms);
    v1=beta*v1;
    locw(rlnfirms) = entry_k;
    srtordr = flipud(sortrows([locw,xothers,oneton],1));
    locw(rlnfirms) = 0;
    locwe = srtordr(:,1); newxe = srtordr(:,2);
    [tempv2,v2] = calcval(locwe,newxe,maxind(srtordr(:,3)==rlnfirms));
    v2 = beta*v2;
    locentr = (v2 - x_entryl - v1) / (x_entryh - x_entryl);
    locentr = min([max([locentr;0]);1]);
  else; locentr = 0;
    end
  if locentr > 0;
    newplace = zeros(rlnfirms,1);
    i=1;
    while i <= rlnfirms;
      newplace(srtordr(i,3)) = i;
      i=i+1;
      end
    end

  v1 = 0; v2 = 0; tempv1 = 0; tempv2 = 0;
  if max(locw) == 0;   % State is all zeros: return value 0

    if locentr < 1;
      [v1, v2] = calcval(locw,newx,1);
      end
    if locentr > 0;
      [tempv1, tempv2] = calcval(locwe,newxe,newplace(1));
      end
    v2 = (1-locentr)*v2 + locentr*tempv2;
    p=0;
    end

  % Get rid of the zeros at the end by finding last non-zero elt.

  i = max([(min(locw) > 0)*rlnfirms;minind(locw)-1]);
  newx = xothers;
  firm = 1;
  while firm <= i;
    if locentr < 1;
      [v1, v2] = calcval(locw,newx,firm);
      end
    if locentr > 0;
      [tempv1, tempv2] = calcval(locwe,newxe,newplace(firm));
      end
    v1 = (1-locentr)*v1 + locentr*tempv1;
    v2 = (1-locentr)*v2 + locentr*tempv2;
    if v1 <= v2; % If value of investing is less, optimal to invest 0

      r = 1.0;
    else; r = 1/(beta*a*(v1-v2));
      end

    % r now contains the value r = (1 - p)^2. => p = 1 - sqrt(r)),
    % where p is the optimal prob. of having k rise, cond. on world

    r = min([max([r;0.00000000001]);1]);
    p = 1 - sqrt(r);
    newx(firm) = p/(a - a * p);
    if locentr > 0;
      newxe(newplace(firm)) = newx(firm);
      end
    firm = firm+1;
    end

  newvalue = profit(ww) - sum(newx) + beta*((1-p)*v2 + p*v1) ...
         - locentr*0.5*((2-locentr)*x_entryl + locentr*x_entryh);

  out1 = newvalue;
  out2 = newx;
  out3 = locentr;


function [out1,out2] = calcval(w,x,i)
% This procedure calculates val = EEEV(.,.,.,.)p(.)p(.)p(.), where E
% represents sums, and this is the calculation of the 4-firm problem
% Vars: w = the vector of everyone's omega
%   x = everyone's investments (rlnfirms of them)
%   i = firm to exclude; do calculation for everyone
%   else's omegas.
% Implicit parameter: oldvalue
% local k_pl1,  % w(i)+1, unless it is at its highest value
%   ii,valA,valB,d,e,probmask,z1,z2,locmask,
%   p_up;  % p_down, p of going up/down for all other firms
% Expand mask to allow for the non-inclusion of the ith plant

  global rlnfirms mask two_n kmax a delta oldvalue
  if rlnfirms > 1;
    if i == 1;
      locmask = [zeros(1,two_n);mask];
    elseif i == rlnfirms;
      locmask = [mask;zeros(1,two_n)];
    else; locmask = [mask(1:i-1,:);zeros(1,two_n);mask(i:rlnfirms-1,:)];
      end
  else; locmask = zeros(1,1);
    end
  x(i) = 0;

  z1 = zeros(rlnfirms,1);
  z2 = kmax*ones(rlnfirms,1);
  p_up = (a .* x) ./ (1 + a .* x);
  % p_down = 1 - p_up;

  valA=0; valB=0;
  ii=1;
  while ii <= two_n;
    % probmask = prodc(locmask(:,ii) .* p_up + (1 - locmask(:,ii)) .* p_down);
    probmask = prod(2 .* locmask(:,ii) .* p_up + 1 - locmask(:,ii) - p_up);
    if probmask == 0;
      ii=ii+1;
    else
      d = flipud(sortrows(flipud(w+locmask(:,ii)),1));
      e = d - 1;

      % Check for evaluation of value fn. at -1

      e = max(([e,z1])')';
      % Check for evaluation of value fn. at kmax+1

      d = min(([d,z2])')';
      valB = valB + ((1-delta)*oldvalue(qencode(d)) ...
        + delta*oldvalue(qencode(e)))*probmask;

      w(i)=w(i)+1;
      d = flipud(sortrows(flipud(w+locmask(:,ii)),1));
      e = d - 1;
      w(i)=w(i)-1;

      % Check for evaluation of value fn. at -1
      e = max(([e,z1])')';
      % Check for evaluation of value fn. at kmax+1
      d = min(([d,z2])')';

      valA = valA + ((1-delta)*oldvalue(qencode(d)) ...
        + delta*oldvalue(qencode(e)))*probmask;
      ii=ii+1;
      end
    end

  out1 = valA;
  out2 = valB;


function [out1] = qencode(ntuple)
% This procedure does a quick encode of any n-tuple given in weakly
% descending order. Encoding is done using a table lookup. Each
% column of the table consists of an n-tuple; the ith column is the ith
% n-tuple to be decoded. The table is stored in the variable "etable".

  global rlnfirms encfirm multfac1 multfac2 etable1 etable2
  if rlnfirms <= encfirm;
    out1 = etable1(sum(ntuple.*multfac1)+1);
  else;
    out1 = etable1(sum(ntuple.*multfac1)+1) ...
      + etable2(sum(ntuple.*multfac2)+1);
    end



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


function [out1] = qdecode(code)
% This procedure does a quick decode of a previously encoded number into
% a weakly descending n-tuple. Decoding is done using a table lookup. Each
% column of the table consists of an n-tuple; the ith column is the ith
% n-tuple to be decoded. The table is stored in the variable "dtable".

  out1 = dtable(:,code);


function [out1] = decode(code)
% This procedure takes a previously encoded number, and decodes it into
% a weakly descending n-tuple (n = nfirms)
% local ntuple,digit,i;

  global rlnfirms binom
  code = code-1;
  ntuple = zeros(rlnfirms,1);
  i = 1;
  while i <= rlnfirms;
    digit = 0;
    while binom(digit+rlnfirms-i+2,digit+2) <= code;
      digit=digit+1;
      end
    ntuple(i) = digit;
    code = code-binom(digit+rlnfirms-i+1,digit+1);
    i = i+1;
    end

  out1 = ntuple;
