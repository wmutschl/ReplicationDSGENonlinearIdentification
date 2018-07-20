% Matrix trace function
% written by Paul Klein for his solab2 package: http://paulklein.se/newsite/codes/codes.php
% See also Gomme and Klein (2011) - Second-order approximation of dynamic models without the use of tensors
% in Journal of Economic Dynamics & Control 35, p.609
%
% Input: matrix X = [X1;X2;...;Xn]
% Output: Y: matrix trace, i.e. trm(X) =  [trace(X1);trace(X2);...;trace(Xn)]
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function Y = tracem(X)
n = size(X,2);
m = size(X,1)/n;
Y = zeros(m,1);
for i=1:m
   Y(i,1)=trace(X((n*(i-1)+1):i*n,1:n));
end
