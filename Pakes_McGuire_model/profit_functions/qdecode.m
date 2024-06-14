function [out1] = qdecode(code,nfirms)
    % Bertrand
    % This procedure takes a previously encoded number, and decodes it into
    % a weakly descending n-tuple (n = nfirms - 1)

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

    %Now convert to format of starting at -7, and jumping by 3's
    
    ntuple = (ntuple.*3 - 7);
    out1 = ntuple;

end
