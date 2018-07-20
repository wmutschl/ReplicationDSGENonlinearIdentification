% Derivative of X(theta)=A(theta)*B(theta)*C(theta)*D(theta) w.r.t to theta
% See Magnus and Neudecker (1999), p. 175
%
% Inputs: Matrices A,B,C,D, and the corresponding derivatives w.r.t theta.
% Output: Derivative of product of A*B*C*D w.r.t. theta
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [dX] = DerivABCD(A,dA,B,dB,C,dC,D,dD)
nparam = size(dA,2);
%   If one or more matrices are left out, they are set to zero
if nargin == 4
    C=speye(size(B,2)); dC=spalloc(numel(C),nparam,0);
    D=speye(size(C,2)); dD=spalloc(numel(D),nparam,0);
elseif nargin == 6
    D=speye(size(C,2)); dD=spalloc(numel(D),nparam,0);
end

dX1 = kron(transpose(D)*transpose(C)*transpose(B),speye(size(A,1)))*dA;
dX2 = kron(transpose(D)*transpose(C),A)*dB;
dX3 = kron(transpose(D),A*B)*dC;
dX4 = kron(speye(size(D,2)),A*B*C)*dD;
dX= dX1+dX2+dX3+dX4;

    
end