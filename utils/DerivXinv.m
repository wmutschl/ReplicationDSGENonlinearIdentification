% Analytical derivative of Z(theta) = inv(X) w.r.t to theta
% See Magnus and Neudecker (1999), p. 183
%
% Inputs:
%       X: matrix X
%       dX: derivative of X w.r.t theta
%       kronflag: use compact kronecker notation or differentials
%   
% Outputs:
%       dZ: Analytical derivative of the kronecker product of X and Y
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [dZ] = DerivXinv(X,dX,kronflag)
if nargin <3
    kronflag = 1;
end
    [n,~] = size(X);    
    ntheta = size(dX,2);
    if kronflag
        dZ = -kron(transpose(X)\speye(n),X\speye(n))*dX;
    else
        dZ = sparse(zeros(n^2,ntheta));
        Xinv = X\speye(n);
        for i=1:ntheta
            dx = reshape(dX(:,i),size(X));            
            dz = -Xinv*dx*Xinv;
            dZ(:,i) = dz(:);
        end
    end
end