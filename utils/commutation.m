% Returns Magnus and Neudecker's commutation matrix of dimensions n by m
%
% Inputs: Dimension of original n times m matrix X
%
% Calls: vec for vectorizing a matrix
%
% Outputs: k: commutation matrix that solves k*vec(X)=vec(X')
%
% Original author: Thomas P Minka (tpminka@media.mit.edu), April 22, 2013
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function k = commutation(n, m)

if nargin < 2
  m = n(2);
  n = n(1);
end

if 0
  % first method
  i = 1:(n*m);
  a = reshape(i, n, m);
  j = vec(transpose(a));
  k = zeros(n*m,n*m);
  for r = i
    k(r, j(r)) = 1;
  end
else
  % second method
  k = reshape(kron(vec(eye(n)), eye(m)), n*m, n*m);
end


function V = vec(A)
    V = A(:); 
end

end