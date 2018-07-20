% This function plots the priors for the parameters
% Based upon Dynare's plot_priors.m and subfunctions
%
% Modified April 16, 2015 by Willi Mutschler (willi@mutschler.eu)

function plot_priors(param)
figurename = 'Priors';
npar = length(param.identif);
[nbplt,nr,nc,lr,lc,nstar] = pltorg(npar);



for plt = 1:nbplt,
    hplt = figure('Name',figurename);
    nstar0 = min(nstar,npar-(plt-1)*nstar);
    for index=1:nstar0
        names = [];
        i = (plt-1)*nstar + index;
        [x,f,abscissa,dens,binf,bsup] = draw_prior_density(i,param);
        nam = param.identif_names(i);
        subplot(nr,nc,index)
        hh = plot(x,f,'-k','linewidth',2);
        set(hh,'color',[0.7 0.7 0.7]);
        box on
        title(nam,'Interpreter','none')
        drawnow
    end
end



function [nbplt,nr,nc,lr,lc,nstar] = pltorg(number)

% Copyright (C) 2004-2008 Dynare Team
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

nrstar = 3;
ncstar = 3;
nstar  = nrstar*ncstar;
nbplt  = 0;
nr     = 0;
nc     = 0;
lr     = 0;
lc     = 0;
if number == 1
    nbplt = 1;
    nr    = 1;
    nc    = 1;
elseif number == 2
    nbplt = 1;
    nr    = 2;
    nc    = 1;
elseif number == 3
    nbplt = 1;
    nr    = 3;
    nc    = 1;
elseif number == 4
    nbplt = 1;
    nr    = 2;
    nc    = 2;
elseif number == 5
    nbplt = 1;
    nr    = 3;
    nc    = 2;
elseif number == 6
    nbplt = 1;
    nr    = 3;
    nc    = 2;    
elseif number == 7
    nbplt = 1;
    nr    = 3;
    nc    = 3;    
elseif number == 8
    nbplt = 1;
    nr    = 3;
    nc    = 3;
elseif number == 9
    nbplt = 1;
    nr    = 3;
    nc    = 3;
else
    if number/nstar == round(number/nstar)
        nbplt = number/nstar;
        nr    = nrstar;
        nc    = ncstar;
        lr    = nr;
        lc    = nc; 
    else
        nbplt = ceil(number/nstar);
        nr    = nrstar;
        nc    = ncstar;
        reste = number-(nbplt-1)*nstar;
        if reste == 1
            lr    = 1;
            lc    = 1;
        elseif reste == 2
            lr    = 2;
            lc    = 1;
        elseif reste == 3
            lr    = 3;
            lc    = 1;
        elseif reste == 4
            lr    = 2;
            lc    = 2;
        elseif reste == 5
            lr    = 3;
            lc    = 2;
        elseif reste == 6
            lr    = 3;
            lc    = 2;    
        elseif reste == 7
            lr    = 3;
            lc    = 3;    
        elseif reste == 8
            lr    = 3;
            lc    = 3;
        end
    end
end

end %pltorg end

function [x,f,abscissa,dens,binf,bsup] = draw_prior_density(indx,param)
% Computes values of prior densities at many points (before plotting)
%
% INPUTS
%    indx          [integer]    Parameter number.
%    param    [structure]  Describes the prior beliefs.
%    
% OUTPUTS
%    x             [double]     Row vector, subset of 'abscissa' such as the density is less than 10
%    f             [double]     Row vector, subset of 'dens' such as the density is less than 10
%    abscissa      [double]     Row vector, abscissa 
%    dens          [double]     Row vector, density
%    binf:         [double]     Scalar, first element of x
%    bsup:         [double]     Scalar, last element of x

pshape  = param.prior.type(find(param.fix==0)); pshape = pshape{indx};
parlow  = param.bound.low(find(param.fix==0)); parlow = parlow{indx};
parup   = param.bound.up(find(param.fix==0)); parup = parup{indx};
par1    = param.prior.par1(find(param.fix==0)); par1 = par1{indx};
par2    = param.prior.par2(find(param.fix==0));par2 = par2{indx};

truncprior = 1e-3;
steps = 200;

if strcmp(pshape,'BETA')% Beta prior
    a = (1-par1)*par1^2/par2^2 - par1;
    b = a*(1/par1 - 1);
    infbound = betainv(truncprior,a,b);
    supbound = betainv(1-truncprior,a,b);
    abscissa = linspace(infbound,supbound,steps);
    dens = betapdf(abscissa,a,b);
elseif strcmp(pshape,'GAMMA')% Gamma prior
    b = par2^2/par1;
    a = par1/b;
    infbound = gaminv(truncprior,a,b);
    supbound = gaminv(1-truncprior,a,b);
    abscissa = linspace(infbound,supbound,steps);
    dens = gampdf(abscissa,a,b);
elseif strcmp(pshape,'NORMAL')% Gaussian prior
    infbound = norminv(truncprior,par1,par2); 
    supbound = norminv(1-truncprior,par1,par2);
    abscissa = linspace(infbound,supbound,steps);
    dens = normpdf(abscissa,par1,par2);  
elseif strcmp(pshape,'INVGAMMA')% Inverse-gamma
    infbound = 1/sqrt(gaminv(1-10*truncprior, par2/2, 2/par1));
    supbound = 1/sqrt(gaminv(10*truncprior, par2/2, 2/par1));
    abscissa = linspace(infbound,supbound,steps);
    dens = exp(lpdfig1(abscissa,par1,par2));  
elseif strcmp(pshape,'UNIFORM')% Uniform prior
    infbound = par1;
    supbound = par2;
    abscissa = linspace(infbound,supbound,steps);
    dens = ones(1, steps) / (supbound-infbound);
end 

binf = abscissa(1);
bsup = abscissa(end);
x = abscissa;
f = dens;
f(find(x<parlow))=0;
f(find(x>parup))=0;
end %draw_prior_density end

function [ldens,Dldens,D2ldens] = lpdfig1(x,s,nu)
% Evaluates the logged INVERSE-GAMMA-1 PDF at x.
%
% X ~ IG1(s,nu) if X = sqrt(Y) where Y ~ IG2(s,nu) and Y = inv(Z) with Z ~ G(nu/2,2/s) (Gamma distribution) 
%
% See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more details.
%
%
% INPUTS     
%    x     [double]  m*n matrix of locations,
%    s     [double]  m*n matrix or scalar, First INVERSE-GAMMA-1 distribution parameters, 
%    nu    [double]  m*n matrix or scalar, Second INVERSE-GAMMA-1 distribution parameters. 
%
% OUTPUTS
%    ldens [double]  m*n matrix of logged INVERSE-GAMMA-1 densities evaluated at x.
%
% SPECIAL REQUIREMENTS
%    none

ldens = -Inf( size(x) ) ;
idx = find( x>0 ) ;    

if length(s)==1
    ldens(idx) = log(2) - gammaln(.5*nu) - .5*nu*(log(2)-log(s)) - (nu+1)*log(x(idx)) - .5*s./(x(idx).*x(idx)) ;
else
    ldens(idx) = log(2) - gammaln(.5*nu(idx)) - .5*nu(idx).*(log(2)-log(s(idx))) - (nu(idx)+1).*log(x(idx)) - .5*s(idx)./(x(idx).*x(idx)) ;
end

if nargout >1 
    if length(s)==1
        Dldens(idx) = - (nu+1)./(x(idx)) + s./(x(idx).^3) ;
    else
        Dldens(idx) = - (nu(idx)+1)./(x(idx)) + s(idx)./(x(idx).^3) ;
    end
end

if nargout == 3 
    if length(s)==1
        D2ldens(idx) =  (nu+1)./(x(idx).^2) - 3*s(idx)./(x(idx).^4) ;
    else
        D2ldens(idx) =  (nu(idx)+1)./(x(idx).^2) - 3*s(idx)./(x(idx).^4) ;
    end
end

end %lpdfig1 end


end%main function end




