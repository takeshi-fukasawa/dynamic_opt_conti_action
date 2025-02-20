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

end
