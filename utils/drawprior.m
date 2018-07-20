% Function that draws from the prior domain, checks whether solution is
% determinate and parameters are not out of bounds.
%
% Input: 
%       DSGE_Model: structure containing DSGE model independent of names
%       Settings: structure containing Settings for order of approximation
% Calls:
%       numeval, SolveModel
%
% Output:
%       pdraw: parameter vector drawn from prior domain which yields determinate solution
%       retcode: 0 if indeterminate solution, 1 if determinate
%
% Based upon Dynare drawprior
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function [pdraw,retcode] = drawprior(DSGE_Model,Settings)
% Initialize
nparam_estim = length(DSGE_Model.param.estim);
shape = DSGE_Model.param.prior.type;
para1 = cell2mat(DSGE_Model.param.prior.par1);
para2 = cell2mat(DSGE_Model.param.prior.par2);
ub = cell2mat(DSGE_Model.param.bound.up);
lb = cell2mat(DSGE_Model.param.bound.low);
pdraw = zeros(nparam_estim,1);
%% Get indices from shapes
fixed_index = find(DSGE_Model.param.fix == 1);
if isempty(fixed_index)
    fixed_draws = 0;
else
    fixed_draws =1;
    shape(fixed_index) = {'FIXED'};
end
beta_index = find(strcmp(shape,'BETA'));
if isempty(beta_index)
    beta_draws = 0;
else
    beta_draws = 1;
end
gamma_index = find(strcmp(shape,'GAMMA'));
if isempty(gamma_index)
    gamma_draws = 0;
else
    gamma_draws = 1;
end
gaussian_index = find(strcmp(shape,'NORMAL'));
if isempty(gaussian_index)
    gaussian_draws = 0;
else
    gaussian_draws = 1;
end
invgamma_index = find(strcmp(shape,'INVGAMMA'));
if isempty(invgamma_index)
    invgamma_draws = 0;
else
    invgamma_draws = 1;
end
uniform_index = find(strcmp(shape,'UNIFORM'));
if isempty(uniform_index)
    uniform_draws = 0;
else
    uniform_draws = 1;
end

%% Draw from distributions
if fixed_draws
    pdraw(fixed_index) = DSGE_Model.param.estim(fixed_index);
end

if uniform_draws
    pdraw(uniform_index) = para1(uniform_index) + (para2(uniform_index)-para1(uniform_index)).*rand(1,length(uniform_index));   
end

if gaussian_draws
    pdraw(gaussian_index) = randn(1,length(gaussian_index)).*para2(gaussian_index) + para1(gaussian_index);
    out_of_bound = find( (pdraw(gaussian_index)'>ub(gaussian_index)) | (pdraw(gaussian_index)'<lb(gaussian_index)));
    while ~isempty(out_of_bound),
        pdraw(gaussian_index(out_of_bound)) = randn(1,length(gaussian_index(out_of_bound))).*para2(gaussian_index(out_of_bound)) + para1(gaussian_index(out_of_bound));
        out_of_bound = find( (pdraw(gaussian_index)'>ub(gaussian_index)) | (pdraw(gaussian_index)'<lb(gaussian_index)));
    end
end

if gamma_draws
    b = para2.^2./para1;
    a = para1./b;
    pdraw(gamma_index) = gamrnd(a(gamma_index),b(gamma_index));
    out_of_bound = find( (pdraw(gamma_index)'>ub(gamma_index)) | (pdraw(gamma_index)'<lb(gamma_index)));
    while ~isempty(out_of_bound),
        pdraw(gamma_index(out_of_bound)) = gamrnd(a(gamma_index(out_of_bound)),b(gamma_index(out_of_bound)));
        out_of_bound = find( (pdraw(gamma_index)'>ub(gamma_index)) | (pdraw(gamma_index)'<lb(gamma_index)));
    end
    clear a b;
end

if beta_draws
    a = (1-para1).*para1.^2./para2.^2 - para1;
    b = a.*(1./para1 - 1);
    pdraw(beta_index) = betarnd(a(beta_index),b(beta_index));
    out_of_bound = find( (pdraw(beta_index)'>ub(beta_index)) | (pdraw(beta_index)'<lb(beta_index)));
    while ~isempty(out_of_bound),
        pdraw(beta_index(out_of_bound)) = betarnd(a(beta_index(out_of_bound)),b(beta_index(out_of_bound)));
        out_of_bound = find( (pdraw(beta_index)'>ub(beta_index)) | (pdraw(beta_index)'<lb(beta_index)));
    end
    clear a b;
end

if invgamma_draws
    pdraw(invgamma_index) = ...
        sqrt(1./gamrnd(para2(invgamma_index)/2,2./para1(invgamma_index)));
    out_of_bound = find( (pdraw(invgamma_index)'>ub(invgamma_index)) | (pdraw(invgamma_index)'<lb(invgamma_index)));
    while ~isempty(out_of_bound),
        pdraw(invgamma_index(out_of_bound)) = ...
            sqrt(1./gamrnd(para2(invgamma_index(out_of_bound))/2,2./para1(invgamma_index(out_of_bound))));
        out_of_bound = find( (pdraw(invgamma_index)'>ub(invgamma_index)) | (pdraw(invgamma_index)'<lb(invgamma_index)));
    end
end

%% Check determinacy of prior draw
[ngra,~,nhes,~,netatilde,~,~,~,~,~,~,nv,ny,~,~,~,~,~,~,~,~,~] = numeval(pdraw,DSGE_Model,Settings.approx,'Numerical');
[~,~,~,~,~,~,~,~,~,~,~,~,retcode] = SolveModel(ngra,nhes,netatilde*transpose(netatilde),nv,ny,Settings.approx);

end % function end