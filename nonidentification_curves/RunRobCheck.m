load('OptionsRobustness'); % This loads the options set in identification_run.m for the An & Schorfheide model, in particular we only check parameters of Taylor rule
addpath('../utils','../models','../models/AnSchorfheide','-begin');
theta0 = DSGE_Model.param.estim;      % Set local point as specified in the GUI  
[Solut0,Deriv0] = EvaluateSparse(theta0,DSGE_Model,Settings.approx,'Analytical'); % Solve Model at theta0
Settings.speed = 'No Speed'; % Do analytical derivatives for innovations once to initialize script files
[G0,Omega0] = gmatrix(Solut0,Deriv0,DSGE_Model,Ident_Test,Settings); % Compute G matrix and spectrum across frequencies for theta0
Omega{1} = Omega0; % store into structre
Settings.speed = 'Speed'; % Do not compute analytical derivatives for innovations anymore, i.e. evaluate script files
index_par = ~DSGE_Model.param.fix; % Index for parameters of Taylor rule

%% 10 points on nonidentification curve reported in Qu and Tkachenko's paper table 1
thet_dir1 = repmat(theta0',10,1); % initialize points on nonidentification curve at theta0 (in particular all other parameters do not change)
load Ni_curves_data %load curve points given in Qu and Tkachenko's paper
thet_dir1(:,index_par)=theta_dir1([1445:1445:14450]'+1,[6 7 8 11]); %select ten equally spaced points along Direction 1
for j = 2:10
    thetaj = thet_dir1(j,:)';      % Set local point
    [Solutj,Derivj] = EvaluateSparse(thetaj,DSGE_Model,Settings.approx,'Analytical'); % Solve model at point from nonidentification curve
    [Gj,Omegaj]=gmatrix(Solutj,Derivj,DSGE_Model,Ident_Test,Settings); % Compute G and spectrum at point from nonidentification curve
    Omega{j} = Omegaj; % Store spectrum
end
maxdev_a1=zeros(10,9); %blanks
maxdev_r1=maxdev_a1;
maxdev_rr1=maxdev_a1;
maxdev_r2=maxdev_a1;
maxdev_a2=maxdev_a1;
maxdev_rr2=maxdev_a1;
windex=maxdev_a1;
for i=1:10
    ad=abs(Omega0-Omega{i});
    rd=ad./abs(Omega0);
    maxdev_a1(i,:)=max(ad); %maximum absolute deviations
    maxdev_r1(i,:)=max(rd); %maximum relative deviations
    for j=1:9
     windex(i,j)=max(find(ad(:,j)==max(ad(:,j))));
     maxdev_rr1(i,j) = rd(windex(i,j),j); % maximum absolute deviations in relative form
    end %measure 2 converts maximum abs deviations in maxdev_a1 and converts them into relative 
end

format('short')
disp('Results: Deviations of spectra across frequencies (direction 1)');
disp('                  Maximum absolute deviations')
disp('              Spectral density matrix element number')
disp('    (1,1)      (2,1)     (3,1)     (2,2)     (3,2)     (3,3)')
disp(maxdev_a1(:,[1,2,3,5,6,9]))
disp('                                       ')
disp('             Maximum absolute deviations in relative form')
disp('              Spectral density matrix element number')
disp('    (1,1)      (2,1)     (3,1)     (2,2)     (3,2)     (3,3)')
disp(maxdev_rr1(:,[1,2,3,5,6,9]))
disp('                                       ')
disp('             Maximum relative deviations ')
disp('              Spectral density matrix element number')
disp('    (1,1)      (2,1)     (3,1)     (2,2)     (3,2)     (3,3)')
disp(maxdev_r1(:,[1,2,3,5,6,9]))


%% General procedure to get points on nonidentification curve given Euler 
[V,D] = eig(G0); % V matrix of eigenvectors, D diagonal matrix with eigenvalues
[~,idxEV] = min(abs(real(diag(D)))); % find index for smallest eigenvalue
c_thet = real(V(:,idxEV)); % eigenvector for smallest eigenvalue
if c_thet(1) < 0; c_thet = -c_thet; end % restrict first element to be positive
h = 1e-5; % step size for Euler method
npoints = 100; % number of points considered
thetajold = theta0;
for j = 1:npoints
    thetaj = theta0; % initialize points on nonidentification curve at theta0 (in particular all other parameters do not change)
    thetaj(find(index_par)) = thetajold(find(index_par)) + c_thet.*h;
    [Solutj,Derivj] = EvaluateSparse(thetaj,DSGE_Model,Settings.approx,'Analytical'); % Solve model at point from nonidentification curve    
    [Gj,Omegaj]=gmatrix(Solutj,Derivj,DSGE_Model,Ident_Test,Settings); % Compute G and spectrum at point from nonidentification curve
    [V,D] = eig(Gj); % V matrix of eigenvectors, D diagonal matrix with eigenvalues
    [~,idxEV] = min(abs(real(diag(D)))); % find index for smallest eigenvalue
    c_thet = real(V(:,idxEV)); % eigenvector for smallest eigenvalue
    if c_thet(1) < 0; c_thet = -c_thet; end % restrict first element to be positive
    S = sprintf('robcheck/iter_%d',j);
    save(S,'thetaj','Omegaj');
    thetajold = thetaj;
end

clear Omega;
irun = 1; 
for j = [1 10 20 30 40 50 60 70 80 90 100] % Select ten equally spaced points
    S = sprintf('robcheck/iter_%d',j);
    load(S)
    Omega{irun} = Omegaj;
    irun = irun+1;
end
maxdev_a1=zeros(10,9); %blanks
maxdev_r1=maxdev_a1;
maxdev_rr1=maxdev_a1;
maxdev_r2=maxdev_a1;
maxdev_a2=maxdev_a1;
maxdev_rr2=maxdev_a1;
windex=maxdev_a1;

for i=1:10
    ad=abs(Omega0-Omega{i});
    rd=ad./abs(Omega0);
    maxdev_a1(i,:)=max(ad); %maximum absolute deviations
    maxdev_r1(i,:)=max(rd); %maximum relative deviations
    for j=1:9
     windex(i,j)=max(find(ad(:,j)==max(ad(:,j))));
     maxdev_rr1(i,j) = rd(windex(i,j),j); % maximum absolute deviations in relative form
    end %measure 2 converts maximum abs deviations in maxdev_a1 and converts them into relative 
end

format('short')
disp('Results: Table 2: Deviations of spectra across frequencies (direction 1)');
disp('                  Maximum absolute deviations')
disp('              Spectral density matrix element number')
disp('    (1,1)      (2,1)     (3,1)     (2,2)     (3,2)     (3,3)')
disp(maxdev_a1(:,[1,2,3,5,6,9]))
disp('                                       ')
disp('             Maximum absolute deviations in relative form')
disp('              Spectral density matrix element number')
disp('    (1,1)      (2,1)     (3,1)     (2,2)     (3,2)     (3,3)')
disp(maxdev_rr1(:,[1,2,3,5,6,9]))
disp('                                       ')
disp('             Maximum relative deviations ')
disp('              Spectral density matrix element number')
disp('    (1,1)      (2,1)     (3,1)     (2,2)     (3,2)     (3,3)')
disp(maxdev_r1(:,[1,2,3,5,6,9]))

