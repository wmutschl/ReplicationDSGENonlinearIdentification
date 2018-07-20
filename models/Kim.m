% Model by Kim (2003) - Functional equivalence between intertemporal and
% multisectoral investment adjustment costs, in: Journal of Economic
% Dynamics and Control 27.4, pp. 533-549
%
% Inputs: 
%       spec: Specification of Model
% Output: 
%       DSGE: structure containing all information about original model symbolically
%
% Modified February 20, 2015 by Willi Mutschler (willi.mutschler@wiwi.uni-muenster.de)

function DSGE = Kim(spec)

%% Model specification (spec)
%    0: no measurement errors, shock on a is normal
%    1: measurement error on all observables (normal), shock on a is normal
%    2: no measurement errors, shock on a is t-distributed

   
%% Define model variables
    % Define shocks
    syms e_a1 e_a2;
    shocks1 = e_a1; shocks2 = e_a2;
    % Define measurement errors
    syms v_c1 v_c2 v_i1 v_i2;
    if spec == 0
        v_c1=0; v_c2=0; v_i1=0; v_i2=0;
        measur_err1 = []; measur_err2 = [];
    elseif spec == 1
        measur_err1 = [v_c1 v_i1]; measur_err2 = [v_c2 v_i2];
    elseif spec == 2
        v_c1=0; v_c2=0; v_i1=0; v_i2=0;
        measur_err1 = []; measur_err2 = [];
    end
    % Add additional variables for stochastic innovations, vector u in the paper
    DSGE.shocks_1 = [shocks1 measur_err1]; DSGE.shocks_2 = [shocks2 measur_err2];
    % Define exogenous and endogenous states, vector x in the paper
    syms k0 k1 a0 a1;
    DSGE.states_0 = [k0 a0]; DSGE.states_1 = [k1 a1];    
    % Define controls, vector y in the paper
    syms c1 c2  i1 i2 obs_c1 obs_c2 obs_i1 obs_i2;
    DSGE.controls_1 = [c1 i1 obs_c1 obs_i1]; DSGE.controls_2 = [c2 i2 obs_c2 obs_i2];
    % Selection Matrix for Observables, note observables are in the last rows of y
    DSGE.SelectMat = [zeros(2,2) eye(2,2)];
 
%% Equilibrium conditions 
    % Define symbols for parameters
    syms alph betae delta thet phi rho_a dumpy

    % Auxiliary parameter and variable    
    s=betae*delta*alph/(1-betae+delta*betae);
    lam1 = (1-s)^thet/c1^(1+thet)/(1+thet);
    lam2 = (1-s)^thet/c2^(1+thet)/(1+thet);
    % Euler equation
    f1 = -lam1*(1+thet)*(i1/s)^thet*(i1/k1/delta)^phi+betae*lam2*(alph*(1+thet)*a1^(1+thet)*k1^(alph*(1+thet)-1)+(1-delta)*(i2/k1/delta)^phi*(1+thet)*(i2/s)^thet);
    % Budget constraint    
    f2 = ((1-s)*(c1/(1-s))^(1+thet) + s*(i1/s)^(1+thet))^(1/(1+thet))- (a0*k0^alph);    
    % Resource constraint
    f3 = k1 -(delta*(i1/delta)^(1-phi)+(1-delta)*k0^(1-phi))^(1/(1-phi));
    % Evolution of technology    
    f4 = -log(a1) + rho_a*log(a0) + e_a1;    
    % Measurement equations
    f5 = obs_c1 - c1 -v_c1;
    f6 = obs_i1 - i1 -v_i1;
    % Auxiliary equation for shock
    f7 = e_a2;
    f8 = v_c2;
    f9 = v_i2;
    %Create function f
    if spec == 0
        DSGE.f = [f1;f2;f3;f4;f5;f6;f7];
    elseif spec == 1
        DSGE.f = [f1;f2;f3;f4;f5;f6;f7;f8;f9];
    elseif spec == 2
        DSGE.f = [f1;f2;f3;f4;f5;f6;f7];    
    end
    
    %% Construct var/cov-matrices of the shocks and measurement errors, sig is the perturbation parameter
    syms sigma_a sigma_c sigma_i dfstudt;
    
    % Make perturbation parameter dependent on std.deviation of one of the shocks
    DSGE.sig=sigma_a;
    
    if spec == 0
        DSGE.eta = sigma_a./DSGE.sig;        
    elseif spec == 1
        DSGE.eta = [sigma_a 0 0;0 sigma_c 0;0 0 sigma_i]./DSGE.sig;
    elseif spec == 2
        DSGE.eta = sqrt(dfstudt/(dfstudt-2))*sigma_a./DSGE.sig;        
    end
        
    DSGE.Sigma = eval(DSGE.sig^2.*DSGE.eta*transpose(DSGE.eta));
    DSGE.etatilde = [zeros(length(DSGE.states_0),length(DSGE.shocks_1)); DSGE.eta];
    
    if spec <2
        DSGE.distribution = 'Gaussian';
        DSGE.dfstudt = [];
    elseif spec ==2
        DSGE.distribution = 'Student-t';
        DSGE.dfstudt = dfstudt;
    end
