% Model by An and Schorfheide (2007) - Bayesian Analysis of DSGE models, in: Econometric Reviews 26(2-4):113-172.
%
% Inputs: 
%       spec: Specification of Model
% Output: 
%       DSGE: structure containing all information about original model symbolically
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)

function DSGE = AnSchorfheide(spec)

%% Model specification (spec)
% 0: no measurement errors
% 1: measurement error on all observables
% 2: no measurement errors, all shocks are t-distributed

%% Define model variables
    % Define shocks, vector u in the paper
    syms e_R1 e_R2 e_g1 e_g2 e_z1 e_z2;
    DSGE.shocks_1 = [e_R1 e_g1 e_z1]; 
    DSGE.shocks_2 = [e_R2 e_g2 e_z2];
    % Define measurement errors
    syms v_YGR1 v_YGR2 v_INFL1 v_INFL2 v_INT1 v_INT2;
    if spec == 0 || spec == 2
        v_YGR1=0; v_INFL1=0; v_INT1=0; v_YGR2=0; v_INFL2=0; v_INT2=0;        
    elseif spec == 1
        % Add measurement errors to shocks
        DSGE.shocks_1 = [DSGE.shocks_1 v_YGR1 v_INFL1 v_INT1];  
        DSGE.shocks_2 = [DSGE.shocks_2 v_YGR2 v_INFL2 v_INT2];
    end
    % Define exogenous and endogenous states, vector x1 and x2 in the paper
    syms y0 y1 R0 R1 g0 g1 z0 z1;
    DSGE.states_0 = [y0 R0 g0 z0]; 
    DSGE.states_1 = [y1 R1 g1 z1];
    % Define controls, vector y in the paper
    syms c1 c2 dy1 dy2 p1 p2 YGR1 YGR2 INFL1 INFL2 INT1 INT2;
    DSGE.controls_1 = [c1 dy1 p1 YGR1 INFL1 INT1]; 
    DSGE.controls_2 = [c2 dy2 p2 YGR2 INFL2 INT2];
    % Selection-Matrix for Observables, note observables are in the last rows of y
    DSGE.SelectMat = [zeros(3,3) eye(3,3)];

%% Equilibrium conditions
    % Define parameter symbolically
    syms tau phi psi1 psi2 rhoR rhog rhoz rA pA gamQ nu cyst;
    syms sigR sigg sigz sig_YGR sig_INFL sig_INT dumpy ;
    % Auxiliary parameters and variables
    pist = exp(pA/400); % page 122, however, AS(2007) code is different
	rrst2 = exp(rA/400);
	bet   = 1/rrst2; % page 122, however, AS(2007) code is different
    
    %phi  = tau*(1-nu)/nu/kap/(pist^2); % eq. (32)
    gst = 1/cyst; %eq. (20)
    
    % Euler equation, eq. (21)
    eq21 = 1 - exp(-tau*c2+tau*c1+R1-rhoz*z1-p2);
    % Phillips curve, eq. (22)
    eq22 = (1-nu)/nu/phi/(pist^2)*(exp(tau*c1)-1)...
        - (exp(p1)-1)*((1-1/2/nu)*exp(p1)+1/2/nu)...
        + bet*(exp(p2)-1)*exp(-tau*c2+tau*c1+dy2+p2);
    % Equilibrium condition, eq. (23)
    eq23 = exp(c1-y1) - exp(-g1) + phi*pist^2*gst/2*(exp(p1)-1)^2;
    % Taylor Rule, eq. (24)
    eq24 = R1 - rhoR*R0 - (1-rhoR)*psi1*p1 - (1-rhoR)*psi2*(y1-g1)-e_R1;
    % Fiscal rule, eq. (25)
    eq25 = g1 - rhog*g0 -e_g1;
    % Evolution of technology, eq (26)
    eq26 = z1 - rhoz*z0 -e_z1;
    % Measurement equations, eq. (38)
    eq38a = YGR1 - gamQ -100*(dy1+z1) - v_YGR1;
    eq38b = INFL1 - pA - 400*p1 - v_INFL1;
    eq38c = INT1 -pA - rA - 4*gamQ - 400*R1 -v_INT1;
    % Output Growth
    eqOutputGrowth  = dy1 - y1 + y0;
    % Auxiliary equations for shocks
    eqshock1 = e_R2;
    eqshock2 = e_g2;
    eqshock3 = e_z2;
    eqmeas1 = v_YGR2;
    eqmeas2 = v_INFL2;
    eqmeas3 = v_INT2;
    %Create function f
    if spec == 0 || spec == 2
        DSGE.f = [eq21;eq22;eq23;eq24;eq25;eq26;eq38a;eq38b;eq38c;...
                eqOutputGrowth;eqshock1;eqshock2;eqshock3];    
    elseif spec == 1
        DSGE.f = [eq21;eq22;eq23;eq24;eq25;eq26;eq38a;eq38b;eq38c;...
                eqOutputGrowth;eqshock1;eqshock2;eqshock3;eqmeas1;eqmeas2;eqmeas3];
    end
%% Construct var/cov-matrices of shocks and measurement errors, sig is the perturbation parameter
    syms sigR sigg sigz sig_YGR sig_INFL sig_INT dfstudt
    
    DSGE.sig=sigz;
    if spec == 0
        DSGE.eta = [sigR 0 0;0 sigg 0;0 0 sigz]./DSGE.sig;
    elseif spec == 1
        DSGE.eta = [sigR 0 0 0 0 0;
                         0 sigg 0 0 0 0;
                         0 0 sigz 0 0 0;
                         0 0 0 sig_YGR 0 0;
                         0 0 0 0 sig_INFL 0;
                         0 0 0 0 0 sig_INT]./DSGE.sig;
    elseif spec == 2
        DSGE.eta = sqrt(dfstudt/(dfstudt-2))*[sigR 0 0;0 sigg 0;0 0 sigz]./DSGE.sig;
    end
    DSGE.Sigma = DSGE.sig^2.*DSGE.eta*transpose(DSGE.eta);
    DSGE.etatilde = [zeros(length(DSGE.states_0),length(DSGE.shocks_1)); DSGE.eta];
    
    if spec <2
        DSGE.distribution = 'Gaussian';
        DSGE.dfstudt = [];
    elseif spec ==2
        DSGE.distribution = 'Student-t';
        DSGE.dfstudt = dfstudt;
    end

