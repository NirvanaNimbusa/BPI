% Sherief Reda (Brown University) and Adel Belouchrani (ENP)
% "Blind Identification of Power Sources in Processors", in IEEE/ACM Design, Automation & Test in Europe, 2017.
% sherief_reda@Brown.edu and adel.belouchrani@enp.edu.dz

function [R,T,pprime]=find_R(fname, n_samples, T_amb)

% This function estimates the thermal matrix R during steady state 
% T=RP
%
% Authors:  S.Reda and A.Belouchrani
% Supported by US NAS Grant 2016, Brown
%

SS = dlmread(fname, '\t');

% steady state total power of cores
c = SS(1:n_samples,1);
% steady state temperatures of all cores
T = SS(1:n_samples, 2:end)';
num_cores = size(T, 1);

T=T-T_amb;

% initialize and run NMF on T
phat=repmat(c'/num_cores, [num_cores 1]);
Rhat=eye(num_cores);
opts.maxit = 20000;  % max # of iterations
opts.tol = 1e-4;     % stopping tolerance
opts.X0 = Rhat;
opts.Y0 = phat;
[Rhat, phat, ~] = nmf(T,num_cores,opts);

% solve permutation through reconstruction 
for iB=1:num_cores
    Tprime(:,:,iB) = Rhat(:,iB)*phat(iB,:); 
    Tm             = mean(abs(Tprime(:,:,iB)),2);
    [~,Imax]    = max(Tm); 
    ki(iB)         = Imax;
end  
 
% Fixing permutation
[p1,ind] = sort(ki);
Rprime   = Rhat(:,ind);
pprime   = phat(ind,:);


% Correct scaling
palpha   = pinv(pprime.')*c;
for i=1:num_cores
    pprime(i, :)=pprime(i,:)*palpha(i);
end
R     = Rprime*diag(1./palpha');    
    
