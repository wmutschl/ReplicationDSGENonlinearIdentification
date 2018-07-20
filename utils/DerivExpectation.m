% This function computes 
% (1) the expectation vector of the observables and its analytical derivative
% (2) the covariance of states and its analytical derivative
% If ChoiceDeriv is set to numerical, then the analytical derivatives are set to zero.
%
% Inputs:
%       Solut: structure containing sparse numerical evaluated solution matrices
%       Deriv: structure containing sparse numerical evaluated derivatives of solution matrices
%       numbers: structure containing size of states, controls and shocks
%       Deriv_type: Analytical or numerical derivatives
%       approx: Order of approximation
%
% Calls:
%       DerivABCD: Derivative of A*B*C*D w.r.t param_identif
%       DerivXkronY: Derivative of kron(X,Y) w.r.t. param_identif
%       vec: vectorize Matrix
%
% Outputs:
%       Ed: expectation vector of observables
%       DEd_Dparam: analytical derivative of expectation vector of observables w.r.t. param_identif (up to order 2)
%       E_xf_xf: expectation of reshape(vec(Sigma_x),nx,nx)
%       DE_xf_xf_dparam: derivative of expectation of vec(Sigma_x)
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [Ed,DEd_Dparam,E_xf_xf,DE_xf_xf_Dparam] = DerivExpectation(Solut,Deriv,numbers,Deriv_type,approx)
SelectMat = Solut.SelectMat;
SS = Solut.SS;                      
prun_A = Solut.prun_A; prun_C = Solut.prun_C;
prun_c = Solut.prun_c; prun_d = Solut.prun_d;
nv=numbers.nv; nx=numbers.nx; nd=numbers.nd;
if strcmp(Deriv_type,'Analytical')
    DSS_Dparam = Deriv.DSS_Dparam;
    nparam = size(DSS_Dparam,2);
    Dprun_A_Dparam = Deriv.Dprun_A_Dparam; Dprun_C_Dparam = Deriv.Dprun_C_Dparam;
    Dprun_c_Dparam = Deriv.Dprun_c_Dparam; Dprun_d_Dparam = Deriv.Dprun_d_Dparam;
else
    nparam=numbers.nparam;
end

% Use notation and formulas of Andreasen et al (2014)
% The auxiliary system is z=[xf xs kron(xf,xf)]
% the mean value of the states in the auxiliary system
IminA = speye(2*nx+nx^2)-prun_A;
E_z     = IminA\prun_c;
E_xf_xf = reshape(E_z(2*nx+1:2*nx+nx^2,1),nx,nx);
% Moments of y: Recall yhat = C*z+d
E_yhat        = prun_C*E_z+prun_d;
% The mean of the states and the controls
if approx == 1
    Ed = SelectMat*SS(nv+1:end);
elseif approx == 2
    Ed = SelectMat*(SS(nv+1:end) + E_yhat);
end
        
% Calculate derivative of expectation of observables w.r.t. param_identif
if strcmp(Deriv_type,'Numerical')
    % If Numerical derivatives, set analytical derivatives to zero
    DEd_Dparam=zeros(nd,nparam); DE_xf_xf_Dparam =zeros(nx^2,nparam);
else % Analytical derivative
    DIminA_Dparam = -Dprun_A_Dparam;
    DinvIminA_Dparam = kron(-transpose(IminA)\speye(size(IminA,2)),IminA\speye(size(IminA,1)))*DIminA_Dparam;
    DE_z_Dparam = DerivABCD(IminA\speye(size(IminA,1)),DinvIminA_Dparam,prun_c,Dprun_c_Dparam);
    DE_xf_xf_Dparam = DE_z_Dparam(2*nx+1:2*nx+nx^2,:);
    DE_yhat_Dparam = DerivABCD(prun_C,Dprun_C_Dparam,E_z,DE_z_Dparam) + Dprun_d_Dparam;
    if approx == 1
        DEd_Dparam = SelectMat*DSS_Dparam(nv+1:end,:);
    elseif approx==2
        DEd_Dparam = SelectMat*(DSS_Dparam(nv+1:end,:) + DE_yhat_Dparam);
    end
end
end