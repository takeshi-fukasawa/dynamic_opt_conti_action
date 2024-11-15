% oct2mat
% gauss2oct

% This is a slightly modified version of markov.g program that allows
% for reading parameters from the configuration file.
% Version 18 uses the Allw x nfirms coding method, and allows for starting
% at any number of firms.
% Written by: Gautam Gowrisankaran
% May 16, 1993

% Modified by Takeshi Fukasawa in April 2024

function [newvalue,newx,iter_info,other_vars] = eql_ma(method,no_entry_exit_spec,spec,c)


% #include pmg.h;
% "**** Computing Dynamic Equilibrium ****";
% #include init.h;
% constants not modifiable by user

global kmax x_entryl x_entryh phi entry_k beta delta a
global binom dtable encfirm etable1 etable2 y mask multfac1 multfac2
global nfirms profit two_n wmax prising
global INV_COST QUAD_INV_COST
global geval_total

kmax = c.KMAX;
x_entryl = c.ENTRY_LOW;
x_entryh = c.ENTRY_HIGH;
phi = c.SCRAP_VAL;
entry_k = c.ENTRY_AT;
rlnfirms = c.MAX_FIRMS;
stfirm = c.START_FIRMS;
beta = c.BETA;
delta = c.DELTA;
a = c.INV_MULT;
INV_COST=c.INV_COST;
QUAD_INV_COST=c.QUAD_INV_COST;

tol = 0.1;  % Tolerance for convergence
newvalue = []; newx = []; oldvalue = []; oldx = []; isentry = [];

% Set up binomial coefficients for decoding/encoding of n-tuples

binom = eye(rlnfirms+kmax+1);
binom = [zeros(rlnfirms+kmax+1,1),binom];
i=2;
while i <= rlnfirms+kmax+1;
  binom(i,2:i) = binom(i-1,2:i) + binom(i-1,1:i-1);
  i=i+1;
  end
encfirm = 3;  % Max. number of firms to encode in table

oneton = 1;
nfirms = stfirm; % max. # of firms at each computation stage

if nfirms > 1;
  nfirms = nfirms - 1;
  wmax = binom(nfirms+1+kmax,kmax+2);

  % Read in data
  % This data is: v (value), x (investment), p (probability of state rising),
  %   isentry

  load(['./middle_output/' 'a.' c.PREFIX '_markov' int2str(nfirms) '.mat']);

  oneton = zeros(nfirms,1);
  i=1;
  while i <= nfirms;
    oneton(i) = i;
    i=i+1;
    end

  if nfirms >= encfirm;
    multfac1 = (kmax+1)^(oneton(1:encfirm)-1);
    nfirms = encfirm;

    % Encode all numbers from 1 to kmax^nfirms

    etable1 = zeros((kmax+1)^nfirms,1);
    i=0;
    while i < rows(etable1);
      msk = [];
      j=kmax+1; k=i;
      while j <= rows(etable1);
        msk = [((mod(k, j))*(kmax+1)/j);msk];
        k = k - (mod(k, j));
        j = j*(kmax+1);
        end
      etable1(i+1) = encode(flipud(sortrows(msk(1:nfirms),1)));
      i=i+1;
      end
    nfirms = stfirm-1;
    if nfirms > encfirm;
      multfac1 = [zeros(nfirms-encfirm,1);multfac1];
      end
    end
  nfirms = stfirm;
  end

while nfirms <= rlnfirms;

% Number of different combinations of competitors positions faced by a firm

  wmax = binom(nfirms+1+kmax,kmax+2);

  disp(sprintf('\nFirms: %d   States: %d\nInitialization ...', nfirms, wmax));

  load(['./middle_output/' 'a.' c.PREFIX '_pr' int2str(nfirms) '.mat'])

  two_n = 2^(nfirms-1);
  dtable = [];

  if nfirms > 1;
    oneton = [oneton;nfirms];
    % Build a mask of all binary numbers from 0 to two_n - 1
    mask = zeros(nfirms-1,two_n);
    i=0;
    while i < two_n;
      msk = [];
      j=2; k=i;
      while j <= two_n;
        if mod(k, j) == 0;
          msk = [0;msk];
        else; k = k - j/2; msk = [1;msk];
          end
        j = j*2;
        end
      mask(:,i+1) = msk(1:nfirms-1);
      i=i+1;
      end
    %print "Mask is " mask;
    end

  % Make a table for quick decoding

  dtable = zeros(nfirms,wmax);
  i=1;
  while i <= wmax;
    dtable(:,i) = decode(i);
    i=i+1;
    end

  % Make a table for quick encoding
  % Fill in multfac1, multfac2, for quick encoding

  if nfirms <= encfirm;
    multfac1 = (kmax+1).^(oneton-1);
    % Encode all numbers from 1 to kmax^nfirms
    etable1 = zeros((kmax+1)^nfirms,1);
    i=0;
    while i < rows(etable1);
      msk = [];
      j=kmax+1; k=i;
      while j <= rows(etable1);
        if isempty(msk) msk = [mod(k, j)*(kmax+1)/j];
        else msk = [mod(k, j)*(kmax+1)/j;msk];
          end
        k = k - (mod(k, j));
        j = j*(kmax+1);
        end
      etable1(i+1) = encode(flipud(sortrows(msk(1:nfirms),1)));
      i=i+1;
      end
  else;
    multfac1 = [0;multfac1];
    multfac2 = [(kmax+1)^(oneton(1:nfirms-encfirm)-1);zeros(encfirm,1)];
    % "Multfac1 is " multfac1;
    % "Multfac2 is " multfac2;
    etable2 = zeros((kmax+1)^(nfirms-encfirm),1);
    i=0;
    while i < rows(etable2);
      msk = [];
      j=kmax+1; k=i;
      while j <= rows(etable2);
        msk = [((mod(k, j))*(kmax+1)/j);msk];
        k = k - (mod(k, j));
        j = j*(kmax+1);
        end
      etable2(i+1) = encode(flipud(sortrows(([msk(1:nfirms-encfirm); ...
        zeros(encfirm,1)]),1)))-1;
      i=i+1;
      end
    end

  % Update values, or define starting values.

  [oldvalue,oldx]=update(newvalue,newx);

  isentry = zeros(wmax,1);
  newx = zeros(wmax,nfirms);
  newvalue = zeros(wmax,nfirms);

  disp(sprintf('Contraction ...'));
  ix = 1;

  %%%%%%%%%%%%%%%%%%%%%%
geval_total=0;
  if 1==1
      [output_spectral,other_vars,iter_info]=...
        spectral_func(@contract,spec,{oldvalue,oldx},...
        method,no_entry_exit_spec);

      newvalue=output_spectral{1};
      newx=output_spectral{2};
      oldvalue=newvalue;
      oldx=newx;

      iter_info.geval_total=geval_total;

      iter_info.feval
  else % Original code (without using the spectral algorithm)
  norm = tol + 1;
  avgnorm = norm;
  while (norm > tol) & (avgnorm > 0.0001*tol);
    [out,other_vars]=contract(oldvalue,oldx,method,no_entry_exit_spec);
    newvalue=out{1};
    newx=out{2};
    norm = max(max(abs(oldvalue - newvalue)));
    avgnorm = mean(mean(abs(oldvalue-newvalue)));

    disp(sprintf('  %2d    Sup norm: %8.4f      Mean norm: %8.4f', ...
      ix, norm, avgnorm));
    ix = ix+1;

    aaaaa = abs(oldvalue-newvalue);
    %normind1 = maxind(max(aaaaa'));
    %normind2 = maxind(max(aaaaa));
    %normcode = (qdecode(normind1))';

    % "Max. elt is: " normind2 "," normcode "; Old value: "
    % oldvalue(normind1,normind2) "; New value: "
    % newvalue(normind1,normind2) "";
  
    oldx = newx; oldvalue = newvalue;
    end

  end%%%%

  % d2 = date;
  % Now find if there is any investment at the highest level.

  w=kmax;
  if nfirms > 1;
    w = [w;zeros(nfirms-1,1)];
    end
  if max(newx(qencode(w,encfirm,etable1,etable2,multfac1,multfac2,nfirms):wmax,1)) > 0;
    disp('Warning: Positive investment recorded at highest efficiency level.')
    disp('Please consider increasing the maximum efficiency level (kmax).')
    end

% Store data in file for inspection
% Store data in file for comparative statics program to read

  prising = a.*newx./(1+a.*newx);
  %save(['a.' c.PREFIX '_markov' int2str(nfirms) '.mat'], ...
  %  'newvalue', 'newx', 'prising', 'isentry')

% disp(sprintf('\n'))
% disp('Value Function (wmax x nfirms)')
% disp([dtable' newvalue])
% disp('Investment (wmax x nfirms)')
% disp([dtable' newx])
% disp('Probability of p rising (wmax x nfirms)'),
% disp([dtable' prising])
% disp('Probability of entry (wmax x nfirms)')
% disp([dtable' isentry])

  nfirms = nfirms+1;
  end

c.EQL_DONE = 1;



function [out,other_vars] = contract(oldvalue,oldx,method,no_entry_exit_spec)
  % This procedure does one iterative step on investment and the value fn
  % local w;
  % First: check for which values of w_s would a firm want to enter

  global wmax isentry
  global a delta kmax mask nfirms two_n
  global encfirm etable1 etable2 multfac1 multfac2
  global beta entry_k x_entryl x_entryh

  isentry=chkentry(oldvalue,oldx,a,delta,kmax,mask,nfirms,two_n,...
    encfirm,etable1,etable2,multfac1,multfac2,beta,entry_k,wmax,x_entryl,x_entryh);

  if no_entry_exit_spec==1%%%%%%
    isentry=isentry*0+1;%entry prob = 1


  end

  % Above is vector of whether firms want to enter, for any w1,...,wn-1

  w = 1;
  while w <= wmax;
    [newx(w,:), newvalue(w,:),diff_mat(w,:)] = optimize(w,oldvalue,oldx,isentry,method,no_entry_exit_spec);
    w=w+1;
  end

  out{1}=newvalue;
  out{2}=newx;
  other_vars=[];
  other_vars.diff_mat=diff_mat;


function [oldvalue,oldx] = update(newvalue,newx)
% This procedure takes the solved newx, newvalue matrix for the nfirms - 1
% problem, and puts them into the nfirms matrices oldx, oldvalue, for use
% as starting values
% local w,i,n,tuple;

  global nfirms wmax
  oldx = zeros(wmax,nfirms);
  oldvalue = zeros(wmax,nfirms);
  if nfirms == 1;
    i=1;
    while i <= wmax;
      oldvalue(i,:) = 1 + 0.1*i;
      i=i+1;
      end
  else;
    w=1;
    while w <= wmax;
      tuple = qdecode(w);
      nfirms = nfirms - 1;
      n = encode(tuple(1:nfirms));
      oldx(w,1:nfirms) = newx(n,1:nfirms);
      oldvalue(w,1:nfirms) = newvalue(n,1:nfirms);
      nfirms = nfirms + 1;
      tuple(nfirms-1) = tuple(nfirms);
      tuple(nfirms) = 0;
      oldvalue(w,nfirms) = oldvalue(encode(tuple),nfirms-1);
      oldx(w,nfirms) = oldx(encode(tuple),nfirms-1);
      w=w+1;
      end
    end



function [out1,out2,diff_temp] = optimize(w,oldvalue,oldx,isentry,method,no_entry_exit_spec)
% This procedure calculates optimal investment, and value fn., for a
% given industry structure w. Thus, a vector nfirms long of each is returned.
% local locw,locwx,locwe,  % Decoded copies of other's omegas w and w/o entry
%   oval,ox, % Old local values
%   entered, % Indicates the probability of an entrant
%   v1,v2,  % v1: value of investing; v2: value of not investing
%   i,j,p,r,tempv1,tempv2,temp, nval,nx; % Returned values of investment, value fn.

  global a beta entry_k nfirms phi profit INV_COST QUAD_INV_COST
  global diff diff_temp entered 
  global encfirm etable1 etable2 multfac1 multfac2
  global delta kmax mask two_n lambda_param
     global geval_total

  locw = qdecode(w);
  locwx = locw;
  oval = oldvalue(w,:)';
  ox = oldx(w,:)';
  nval = zeros(nfirms,1);
  nx = zeros(nfirms,1);
  diff=zeros(nfirms,1);
  diff_temp=zeros(nfirms,1);

  if no_entry_exit_spec==0%%%%%%%%%%%
    % Find out which firms want to exit

    i = (min(oval) == phi)*(minind(oval)-1) + (min(oval) > phi)*nfirms;

    % Replace efficiency levels of exitors with zero

    if i < nfirms;
      locwx(i+1:nfirms) = zeros(nfirms-i,1);
    end
  end


  % Figure out the probability of entry

  entered = isentry(qencode(flipud(sortrows(flipud(locwx),1)),encfirm,etable1,etable2,multfac1,multfac2,nfirms));

  locwe = locwx;
  locwe(nfirms) = entry_k;

  % Now calculate the optimal policies for this industry structure, given that
  % entry and exit are as specified.

  j=1;
  while j <= nfirms;
    if locw(j) == 0;
      nval(j:nfirms) = phi*ones(nfirms-j+1,1);
      break;
    end

    v1=0; v2=0;
    if entered < 1;

      % First: Calculate v, without entry

      [v1, v2] = calcval(j,locwx,ox,locw(j),oldvalue,...
        a,delta,kmax,mask,nfirms,two_n,...
          encfirm,etable1,etable2,multfac1,multfac2);
    end

    if entered > 0;

      % A firm wants to enter with positive probability

      [tempv1, tempv2] = calcval(j,locwe,ox,locw(j),oldvalue,...
        a,delta,kmax,mask,nfirms,two_n,...
          encfirm,etable1,etable2,multfac1,multfac2);
      v1 = entered*tempv1 + (1-entered)*v1;
      v2 = entered*tempv2 + (1-entered)*v2;
    end


    if method=="PM" & QUAD_INV_COST==0
        % Calculate values for firm, given that it is not leaving

        if v1 <= v2; % Avoid division by zeros
            r = 1.0;
        else; r = 1.0/(beta*a*(v1-v2));
        end

        % r now contains the value r = (1 - p)^2. => p = 1 - sqrt(r)),
        % where p is the optimal prob. of having k rise, cond. on world

        r = min([max([r;0.0000000000001]);1]);
        p = 1.0 - sqrt(r); % based on the optimality condition

        nx(j) = p/(a - a * p);
    
        geval_total=geval_total+1;


     elseif method=="PM" & QUAD_INV_COST>0

        %%% Nonlinear optimization
        %options = optimset('Display','off');
        %x_min=0;x_max=1;%%%
        %x_sol=fminbnd(@Q_func,x_min,x_max,options,v1,v2,a,beta,INV_COST,QUAD_INV_COST);
        
        [x_sol,exitflag,n_iter,geval]=csolve('FOC_func',ox(j),[],0.000001,5,...
              v1,v2,a,beta,INV_COST,QUAD_INV_COST); 
        geval_total=geval_total+geval;

        %%% Use fmincon or fminunc => Slower...
        %x_j_init=ox(j);
        %options=[];
        %A=[];b=[];Aeq=[];beq=[];
        %lb=0;ub=[];nonlcon=[]; 
        %options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',true);
        %x_sol=fmincon(@Q_func,x_j_init,A,b,Aeq,beq,lb,ub,nonlcon,options,v1,v2,a,beta,INV_COST,QUAD_INV_COST);
        
        %options = optimoptions('fminunc','Display','off');
        %x_sol=fminunc(@Q_func,x_j_init,options,v1,v2,a,beta,INV_COST,QUAD_INV_COST);
        
        if x_sol<1e-4
            x_sol=0;% avoid warning message concerning the range of states
        end

        nx(j)=x_sol;

        p=(a.*ox(j))./(1+a.*ox(j));

    elseif method=="gradient"

        [diff(j),p]=FOC_func(ox(j),v1,v2,a,beta,INV_COST,QUAD_INV_COST);
        if abs(ox(j))<=1e-16 & diff(j)<0
            diff(j)=0;
        end

        nx(j)=ox(j)+lambda_param*diff(j);

        geval_total=geval_total+1;

    end
    
    % Now calculate the value from staying in
    % Ask: given this optimal investment level, will there be exit?


    %%% Compute nval based on updated x (nx)
    %%nval(j) = profit(w,j) - INV_COST.*nx(j)-QUAD_INV_COST.*nx(j).^2 + beta*(v1*p + v2*(1-p));

    %%% Compute nval based on initial x (ox)
    nval(j) = profit(w,j) - INV_COST.*ox(j)-QUAD_INV_COST.*ox(j).^2 + beta*(v1*p + v2*(1-p));

    if nval(j) <= phi & no_entry_exit_spec==0
      nval(j) = phi;
      nx(j) = 0;
    end

    if (j < nfirms) & (nval(j) == phi) & no_entry_exit_spec==0
      nval(j+1:nfirms) = ones(nfirms-j,1) * phi;
      break;
    end

    %%%%ox(j) = nx(j); %% Sequential update%%%%

    if no_entry_exit_spec==0
        locwx(j) = (nval(j) > phi)*locw(j);
    else
        locwx(j) = locw(j);
        locwe(j) = locwx(j);
        j=j+1;
    end
  end% while


  out1 = nx';
  out2 = nval';


function [isentry] = chkentry(oldvalue,oldx,a,delta,kmax,mask,nfirms,two_n,...
  encfirm,etable1,etable2,multfac1,multfac2,beta,entry_k,wmax,x_entryl,x_entryh)
% This procedure calculates for which value of other people's omegas, would
% a firm want to enter, given that there is room in the market for
% a firm to enter
% local w,locw,v1,vgarbage,
%   val; % Value from entering

  isentry=zeros(wmax,1);
  w = 1;
  while w <= wmax;
    locw = qdecode(w);
    if locw(nfirms) == 0;
      [vgarbage,v1] = calcval(nfirms,locw,oldx(w,:)',entry_k,oldvalue,...
        a,delta,kmax,mask,nfirms,two_n,...
          encfirm,etable1,etable2,multfac1,multfac2);
      val = beta * v1;

      % print val-x_entry;

      isentry(w) = (val - x_entryl) / (x_entryh - x_entryl);
      end
    w=w+1;
    end
  isentry = min([isentry,ones(wmax,1)]')';
  isentry = max([isentry,zeros(wmax,1)]')';



function [out1,out2] = calcval(place,w,x,k,oldvalue,...
    a,delta,kmax,mask,nfirms,two_n,...
    encfirm,etable1,etable2,multfac1,multfac2)
% This procedure calculates val = EEEV(.,.,.,.)p(.)p(.)p(.), where E
% represents sums, and this is the calculation of the 4-firm problem
% Vars: place = place of own omega, for calculating value function (v)
%       w = the vector of omegas; already decoded
%       x = the vector of investments (nfirms of them)
% For efficiency reasons, it outputs the following vector:
% [ calcval(k_v+1,w,x,oldvalue), calcval(k_v,w,x,oldvalue) ]
% local i,valA,valB,d,e,probmask,z1,z2,locmask,
%   p_up,  % p_down, p of going up/down for all other firms
%   temp,
%   pl1,justone;

  z1 = zeros(nfirms,1);
  z2 = kmax*ones(nfirms,1);

  % Expand mask to allow for the non-inclusion of the ith plant

  if nfirms > 1;
    if place == 1;
      locmask = [zeros(1,two_n);mask];
    elseif place == nfirms;
      locmask = [mask;zeros(1,two_n)];
    else; locmask = [mask(1:place-1,:);zeros(1,two_n);mask(place:nfirms-1,:)];
      end
  else; locmask = zeros(1,1);
    end
  x(place) = 0;
  w(place) = k;
  justone = zeros(nfirms,1);
  justone(place) = 1;
  p_up = (a .* x) ./ (1 + a .* x);
  % p_down = 1 - p_up;
  valA=0; valB=0;
  i=1;

  while i <= two_n;
    % probmask = prod(mask(:,i) .* p_up + (1 - mask(:,i)) .* p_down);
    probmask = prod(2 .* locmask(:,i) .* p_up + 1 - locmask(:,i) - p_up);
    d = w+locmask(:,i);
    temp = flipud(sortrows(flipud([d,justone]),1));
    d = temp(:,1);
    e = d - 1;

    % Check for evaluation of value fn. at -1
    e = max(([e,z1])')';
    % Check for evaluation of value fn. at kmax+1
    d = min(([d,z2])')';
    pl1 = maxind(temp(:,2));% sum(d(1:place)>=k) + sum(d(place:nfirms)>k);

    valB = valB + ((1-delta)*oldvalue(qencode(d,encfirm,etable1,etable2,multfac1,multfac2,nfirms),pl1) ...
            + delta*oldvalue(qencode(e,encfirm,etable1,etable2,multfac1,multfac2,nfirms),pl1))*probmask;

    d = w+locmask(:,i)+justone;
    temp = flipud(sortrows(flipud([d,justone]),1));
    d = temp(:,1);
    e = d - 1;

    % Check for evaluation of value fn. at -1
    e = max(([e,z1])')';
    % Check for evaluation of value fn. at kmax+1
    d = min(([d,z2])')';
    pl1 = maxind(temp(:,2)); %sum(e(1:place)>=k) + sum(e(place:nfirms)>k);

    valA = valA + ((1-delta)*oldvalue(qencode(d,encfirm,etable1,etable2,multfac1,multfac2,nfirms),pl1) ...
            + delta*oldvalue(qencode(e,encfirm,etable1,etable2,multfac1,multfac2,nfirms),pl1))*probmask;
    i=i+1;
    end

  out1 = valA;
  out2 = valB;


function [out1] = encode(ntuple)
% This procedure takes a weakly descending n-tuple (n = nfirms), with
% min. elt. 0, max. elt. kmax, and encodes it into an integer
% local code,digit,i;

  global binom nfirms
  code = 1; % Coding is from 1 to wmax
  i = 1;
  while i <= nfirms;
    digit = ntuple(i);
    code = code + binom(digit+nfirms+1-i,digit+1);
    i=i+1;
    end

  out1 = code;


function [out1] = qencode(ntuple,encfirm,etable1,etable2,multfac1,multfac2,nfirms)
% This procedure does a quick encode of any n-tuple given in weakly
% descending order. Encoding is done using a table lookup. Each
% column of the table consists of an n-tuple; the ith column is the ith
% n-tuple to be decoded. The table is stored in the variable "etable".

  if nfirms <= encfirm;
    out1 = etable1(sum(ntuple.*multfac1)+1);
  else;
    out1 = etable1(sum(ntuple.*multfac1)+1) ...
      + etable2(sum(ntuple.*multfac2)+1);
    end



function [out1] = qdecode(code)
% This procedure does a quick decode of a previously encoded number into
% a weakly descending n-tuple. Decoding is done using a table lookup. Each
% column of the table consists of an n-tuple; the ith column is the ith
% n-tuple to be decoded. The table is stored in the variable "dtable".

  global dtable

  out1 = dtable(:,code);


function [out1] = decode(code)
% This procedure takes a previously encoded number, and decodes it into
% a weakly descending n-tuple (n = nfirms)
% local ntuple,digit,i;

  global binom nfirms
  code = code-1;
  ntuple = zeros(nfirms,1);
  i = 1;
  while i <= nfirms;
    digit = 0;
    while binom(digit+nfirms-i+2,digit+2) <= code;
      digit=digit+1;
      end
    ntuple(i) = digit;
    code = code-binom(digit+nfirms-i+1,digit+1);
    i = i+1;
    end

  out1 = ntuple;

