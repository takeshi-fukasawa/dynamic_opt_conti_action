
function FOC = cfunk(p)
%   used for quality competition profit function
   global egw M mc w sigma

   egw = eg(w);
   n = egw.*exp(-p);
   denom = 1.0 + sum(n);
   sigma = n./denom;
  FOC=-(p-mc).*sigma.*(1-sigma) + sigma;
  %FOC=-(p-mc).*(1-sigma) + 1;
end