% Solves for the first and second-order accurate approximation of the solution to the DSGE model
% For details see Gomme and Klein (2011): Second-order approximations of dynamic models without the use of tensors
% These things are edited from original solab2.m:
    % - Notation is adapted to the paper.
    % - More inputs: approx determines the order of approximation.
    % - More outputs: auxiliary matrices are added to the output
    % - kx and ky are multiplied by 2, so they are equivalent to Schmitt-Groh� and Uribe's hss and gss.
    % - Calculations of second order solution matrices are slightly different to be similar to the paper
    % - solab.m is embedded
%          
% Inputs: 
%       gra: an (nv+ny) by 2*(nv+ny) matrix (the gradient of f), evaluated at the steady state
%       hes: an 2*(nv+ny)^2 by 2*(nv+ny) matrix (the Magnus-Neudecker-Hessian of f), evaluated at the steady-state
%       eta_etaT: an nv by nv matrix: the variance of shocks and measurement errors in transition equation (without perturbation parameter)
%       nv=nx+nu: length of state vector and auxiliary vector for shocks
%       ny: length of control vector
%       approx: order of approximation
%
% Calls: Embedded function solab.m, tracem: Matrix trace
%
% Outputs:
%       gv: numerical evaluated first-order solution matrix for controls
%       gvv: numerical evaluated second-order solution matrix for controls
%       gSS: numerical evaluated second-order solution matrix for controls
%       hv: numerical evaluated first-order solution matrix for states
%       hvv: numerical evaluated second-order solution matrix for states
%       hSS: numerical evaluated second-order solution matrix for states
%       solM: [hx; gx*hx; eye(nx); gx]
%       solN: ([eye(nx); gx; zeros(n,nx)]
%       solQ: [kron(hx',kron(f2,hx')) + kron(eye(nx),kron(f4,eye(nx))), kron(eye(nx),(kron(f1,eye(nx))+kron(f2*gx,eye(nx))))]
%       retcode: 0 if indeterminate solution, 1 if determinate
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [gv,gvv,gSS,hv,hvv,hSS,solM,solN,solQ,solR,solS,solU,retcode] = SolveModel(gra,hes,eta_etaT,nv,ny,approx)

n = nv + ny;   % Number of variables

%% Approximation to first order
A = -gra(:,1:n); % A = -[f1 f2]
B = gra(:,n+1:end); % B = [f3 f4]
[gv,hv] = solab(A,B,nv); % First-order solution matrix using Generalized Schur

if isempty(gv)
    % indeterminate solution
    retcode = 0; 
    gvv=[];hvv=[];hSS=[];gSS=[];solM=[];solN=[];solQ=[];solR=[];solS=[];solU=[];
else
    retcode = 1; % determinate solution
    if approx == 1 % Set all other outputs to zero if linear approximation to first order
        gvv=zeros(ny*nv,nv);
        hvv=zeros(nv*nv,nv);
        hSS=zeros(nv,1);
        gSS=zeros(ny,1);
        solM=zeros(2*n,nv);
        solN=zeros(2*n,nv);
        solQ=zeros(n*nv^2,n*nv^2);
        solR=zeros(n*nv,nv);
        solS=zeros(n,n);
        solU=zeros(n,1);
        
    %% Approximation to second order
    elseif approx == 2 
        f1 = gra(:,1:nv);
        f2 = gra(:,nv+1:n);
        f4 = gra(:,n+nv+1:2*n);

        solM = [hv;gv*hv;eye(nv);gv];
        solR = kron(eye(n),solM')*hes*solM;
        solQ1 = kron(hv',kron(f2,hv')) + kron(eye(nv),kron(f4,eye(nv)));
        solQ2 = kron(eye(nv),(kron(f1,eye(nv))+kron(f2*gv,eye(nv))));
        solQ = [solQ1 solQ2];
        vec_gXX_hXX = -solQ\solR(:);
        gvv = reshape(vec_gXX_hXX(1:nv^2*ny),nv*ny,nv);
        hvv = reshape(vec_gXX_hXX(nv^2*ny+1:end),nv^2,nv);

        solS = [f1+f2*gv f2+f4];
        solN = [eye(nv);gv;zeros(n,nv)];
        solU = f2*tracem(kron(eye(ny),eta_etaT)*gvv)+tracem(kron(eye(n),solN')*hes*solN*eta_etaT);
        hss_gss = -solS\solU;
        hSS = hss_gss(1:nv);
        gSS = hss_gss(nv+1:end);
    end
end

% Embedded Solab function
% Solves for the recursive representation of the stable solution to a system of linear difference equations.
% Written by Paul Klein, for details see Gomme and Klein (2011): Second-order approximations of dynamic models without the use of tensors
% See also An & Schorfheide (2007)'s code to control for warning statements
function [f,p] = solab(a,b,nk)
realsmall=1e-7; 
    try
        [s,t,q,z] = qz(a,b);                % upper triangular factorization of the matrix pencil b-za
    catch
        warning('Solab:noexist','Error using qz');
	f = [];
	p = [];
	return;
    end    
zxz = sum((abs(diag(s))<realsmall) & (abs(diag(t))<realsmall));
if ~(~zxz)
	warning('Solab:noexist','Coincident zeros');
	f = [];
	p = [];
	return;
end

    [s,t,~,z] = ordqz(s,t,q,z,'udo');   % reordering of generalized eigenvalues with the block inside the unit circle in the upper left
if abs(t(nk,nk))>abs(s(nk,nk))
	warning('Solab:noexist','No equilibrium exists.');
	f = [];
	p = [];
	return;
elseif abs(t(nk+1,nk+1))<abs(s(nk+1,nk+1))
	warning('Solab:indeterminacy','Indeterminacy.');
	f = [];
	p = [];
	return;
else
	lastwarn('');
end
    
    z21 = z(nk+1:end,1:nk);
    z11 = z(1:nk,1:nk);

    if rank(z11)<nk
        warning('Solab:noexist','Invertibility condition violated')
        f = [];
        p = [];
        return;
    end

    z11i = z11\eye(nk);
    s11 = s(1:nk,1:nk);
    t11 = t(1:nk,1:nk);

%     if abs(t(nk,nk))>abs(s(nk,nk)) | abs(t(nk+1,nk+1))<abs(s(nk+1,nk+1));
%        warning('Wrong number of stable eigenvalues.');
%     end

    dyn = s11\t11;

    f = real(z21*z11i); % The real function takes away very small imaginary parts of the solution
    p = real(z11*dyn*z11i);
end

end %function end