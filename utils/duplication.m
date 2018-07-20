% Duplication Matrix as defined by Magnus and Neudecker (2002), p.49
%
% Inputs:
%       n: size of vector
% Outputs:
%       d: Duplication matrix
%       QPinv: Moore-Penrose inverse of QP
%
% Author: Thomas P Minka (tpminka@media.mit.edu)
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function d = duplication(n)

if 0
  % first method
  a = zeros(n);
  k = 1;
  for j = 1:n
    for i = 1:n
      if i >= j
	a(i,j) = k;
	k = k + 1;
      else
	a(i,j) = a(j,i);
      end
    end
  end
else
  % second method
  a = tril(ones(n));
  i = find(a);
  a(i) = 1:length(i);
  a = a + transpose(tril(a,-1));
end
j = a(:);

m = n*(n+1)/2;
d = zeros(n*n,m);
for r = 1:size(d,1)
  d(r, j(r)) = 1;
end
