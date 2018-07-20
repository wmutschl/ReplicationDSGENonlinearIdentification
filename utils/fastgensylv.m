% Solves the Sylvester equation A * X + B * X * C + D = 0 for X.
% From DYNARE

function X = fastgensylv(A, B, C, D, tol,maxit,X0)

%@info:
%! @deftypefn {Function File} {[@var{X1}, @var{info}] =} fastgensylv (@var{A},@var{B},@var{C},@var{tol},@var{maxit},@var{X0})
%! @anchor{fastgensylv}
%! @sp 1
%! Solves the Sylvester equation A * X + B * X * C + D = 0 for X.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Square matrix of doubles, n*n.
%! @item B
%! Square matrix of doubles, n*n.
%! @item C
%! Square matrix of doubles, n*n.
%! @item tol
%! Scalar double, tolerance parameter.
%! @item maxit
%! Integer scalar, maximum number of iterations.
%! @item X0
%! Square matrix of doubles, n*n, initial condition.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item X
%! Square matrix of doubles, n*n, solution of the matrix equation.
%! @item info
%! Scalar integer, if nonzero the algorithm failed in finding the solution of the matrix equation.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2012 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if size(A,1)~=size(D,1) || size(A,1)~=size(B,1) || size(C,2)~=size(D,2)
    error('fastgensylv:: Dimension error!')
end

if nargin<7 || isempty(X0)
    X = zeros(size(A,2),size(C,1));
elseif nargin==7 && ~isempty(X0)
    X = X0;
end

kk = 0;
cc = 1+tol;

iA = inv(A);
Z = - (B * X * C + D);

while kk<=maxit && cc>tol
    X = iA * Z;
    Z_old = Z;
    Z = - (B * X * C + D);
    cc = max(sum(abs(Z-Z_old)));
    kk = kk + 1;
end

if kk==maxit && cc>tol
    error(['fastgensylv:: Convergence not achieved in fixed point solution of Sylvester equation after ' int2str(maxit) ' iterations']);
end




% function X = fastgensylv(A, B, C, D)
% Solve the Sylvester equation:
% A * X + B * X * C + D = 0
% INPUTS
%   A
%   B
%   C
%   D
%   block : block number (for storage purpose) 
%   tol : convergence criteria
% OUTPUTS
%   X solution
%    
% ALGORITHM
%   fixed point method
%   MARLLINY MONSALVE (2008): "Block linear method for large scale
%   Sylvester equations", Computational & Applied Mathematics, Vol 27, n°1,
%   p47-59
%   ||A^-1||.||B||.||C|| < 1 is a suffisant condition:
%    - to get a unique solution for the Sylvester equation
%    - to get a convergent fixed-point algorithm
%
% SPECIAL REQUIREMENTS
%   none.  
