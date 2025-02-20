function wret = eg(w)
  % used for quality competition profit function
  % Calculates e^g(w)
  global wstar
   n_var=size(w(:),1);
   wret = zeros(n_var,1);
   i=1;
   while i <= n_var
     if w(i) <= wstar;
       wret(i) = exp(w(i));
    else; wret(i) = exp(wstar)*(2.0-exp(-(w(i)-wstar)));
    end
    i=i+1;
  end
end
