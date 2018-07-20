% Analytical derivative of Z(theta) = kron(X(theta),Y(theta)) w.r.t to theta
% See Magnus and Neudecker (1999), p. 185
%
% Inputs:
%       X and Y: matrices X and Y
%       dX and dY: derivatives of X and Y w.r.t theta
%       kronflag: 1 use formula, else use differential
%
% Calls: commutation: get commutation matrix; vec: vectorize matrix
%   
% Outputs:
%       dZ: Analytical derivative of the kronecker product of X and Y
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [dZ] = DerivXkronY(X,dX,Y,dY,kronflag)
if nargin <5
    kronflag = 1;
end
    [n,q] = size(X);
    [p,r] = size(Y);
    ntheta = size(dX,2);
    if kronflag
        dZ = kron(kron(speye(q),sparse(commutation(r,n))),speye(p))*(...
            kron(speye(n*q),Y(:))*dX + kron(X(:),speye(p*r))*dY);
    else
        %dZ = spalloc(n*p*q*r,ntheta,n*p*q*r*ntheta);
        for i=1:ntheta
            dx = reshape(dX(:,i),size(X));
            dy = reshape(dY(:,i),size(Y));
            dz = kron(dx,Y) + kron(X,dy);
            dZ(:,i) = dz(:);            
        end
    end
end