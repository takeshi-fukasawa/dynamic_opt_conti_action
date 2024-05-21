function [x,weight]=gausshermi(n)
%% Gauss-Hermite quadrature 

p=hermipol(n);
%Roots
x=roots(p(n+1,:));

for i=1:n
   weight(i,1)=(2.^(n-1)*(factorial(n)).*sqrt(pi))./(n.^2.*(polyval(p(n,1:n),x(i))).^2);
end
return


function p=hermipol(n)
k=2;
p(1,1)=1;
p(2,1:2)=[2 0];
for k=2:n
   p(k+1,1:k+1)=2*[p(k,1:k) 0]-2*(k-1)*[0 0 p(k-1,1:k-1)];
end
return