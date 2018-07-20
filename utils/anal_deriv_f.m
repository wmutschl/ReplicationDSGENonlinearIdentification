% Computes analytically Jacobian and Magnus-Neudecker-Hessian of the function f(yp,y,xp,x) with respect to x, y, xp, and yp.
% Note: This program requires MATLAB's Symbolic Math Toolbox
% 
% Inputs:
%       f: symbolic equilibrium conditions of DSGE model
%       xp: symbolic vector of states in period t+1
%       yp: symbolic vector of controls in period t+1
%       x: symbolic vector of states in period t
%       y: symbolic vector of controls in period t
%       approx: order of approximation
%           If approx is set at a value different from 2, the program delivers the Jacobian of f and sets the Hessian at zero. 
%           If approx equals 2, the program returns Jacobian and Hessian of f. The default value of approx is 2. 
% Output:
%       gra: symbolic gradient of f
%       hes: symbolic Magnus-Neudecker-Hessian of f
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [gra,hes]=anal_deriv_f(f,xp,yp,x,y,approx) 

if nargin==5
approx=2;
end

% Compute Jacobian of f
gra = jacobian(f,[xp yp x y]);

% Compute Magnus-Neudecker-Hessian of f
if approx==2
    graT = transpose(gra);
    hes = jacobian(graT(:),[xp yp x y]);
else
    hes = zeros(numel(gra),length([xp yp x y]));
end