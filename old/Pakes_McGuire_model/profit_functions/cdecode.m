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
end
