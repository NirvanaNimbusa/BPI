function [X,Y,out] = nmf(M,r,opts)
% <http://www.caam.rice.edu/~optimization/bcu/nmf/nmf.m Download this m-file>
%  nmf: nonnegative matrix factorization by block-coordinate update
%   minimize 0.5*||M - X*Y||_F^2
%   subject to X>=0, Y>=0
%
%  input:
%       M: input nonnegative matrix
%       r: rank (X has r columns and Y has r rows)
%       opts.
%           tol: tolerance for relative change of function value, default: 1e-4
%           maxit: max number of iterations, default: 500
%           maxT: max running time, default: 1e3
%           rw: control the extrapolation weight, default: 1
%           X0, Y0: initial X and Y, default: Gaussian random matrices
%  output:
%       X, Y: solution nonnegative factors
%       out.
%           iter: number of iterations
%           hist_obj: objective value at each iteration
%           relerr: relative changes at each iteration
%
% More information can be found at:
% http://www.caam.rice.edu/~optimization/bcu/

%% Parameters and defaults
if isfield(opts,'tol'),    tol = opts.tol;     else tol = 1e-4;   end
if isfield(opts,'maxit'),  maxit = opts.maxit; else maxit = 500;  end
if isfield(opts,'maxT'),   maxT = opts.maxT;   else maxT = 1e3;   end
if isfield(opts,'rw'),     rw = opts.rw;       else rw = 1;       end

%% Data preprocessing and initialization
[m,n] = size(M);

if isfield(opts,'X0'), X0 = opts.X0; else X0 = max(0,randn(m,r)); end
if isfield(opts,'Y0'), Y0 = opts.Y0; else Y0 = max(0,randn(r,n)); end

Mnrm = norm(M,'fro');

X0 = X0/norm(X0,'fro')*sqrt(Mnrm);
Y0 = Y0/norm(Y0,'fro')*sqrt(Mnrm);
Xm = X0; Ym = Y0; 
Yt = Y0'; Ys = Y0*Yt; MYt = M*Yt;
obj0 = 0.5*Mnrm*Mnrm;

nstall = 0; t0 = 1;
Lx = 1; Ly = 1;

%% Iterations of block-coordinate update 
%
%  iteratively updated variables:
%       Gx, Gy: gradients with respect to X, Y
%       X, Y: new updates
%       Xm, Ym: extrapolations of X, Y
%       Lx, Lx0: current and previous Lipschitz bounds used in X-update
%       Ly, Ly0: current and previous Lipschitz bounds used in Y-update
%       obj, obj0: current and previous objective values
%
%  cached computation:
%       Xs = X'*X, Ys = Y*Y' 

start_time = tic;
fprintf('Iteration:     ');

for k = 1:maxit
    fprintf('\b\b\b\b\b%5i',k);    
    
    % --- X-update ---
    Lx0 = Lx; Lx  = norm(Ys);  % save and update Lipschitz bound for X
    Gx = Xm*Ys - MYt; % gradient at X=Xm
    X = max(0, Xm - Gx/Lx);
    Xt = X'; Xs = Xt*X; % cache useful computation
    
    % --- Y-update ---
    Ly0 = Ly; Ly = norm(Xs);  % save and update Lipschitz bound for Y
    Gy = Xs*Ym - Xt*M; % gradient at Y=Ym
    Y = max(0, Ym - Gy/Ly);
    Yt = Y'; Ys = Y*Yt; MYt = M*Yt;  % cache useful computation
    
    % --- diagnostics, reporting, stopping checks ---
    obj = 0.5*(sum(sum(Xs.*Ys)) - 2*sum(sum(X.*MYt)) + Mnrm*Mnrm);
    
    % reporting
    out.hist_obj(k) = obj;
    out.relerr1(k) = abs(obj-obj0)/(obj0+1);
    out.relerr2(k) = sqrt(2*obj)/Mnrm;
    
    % stall and stopping checks
    crit = (out.relerr1(k)<tol);
    if crit; nstall = nstall+1; else nstall = 0; end
    if nstall >= 3 || out.relerr2(k) < tol, break; end
    if toc(start_time) > maxT; break; end;
    
    % --- correction and extrapolation ---
    t = (1+sqrt(1+4*t0^2))/2;
    if obj>=obj0 
        % restore to previous X,Y, and cached quantities for nonincreasing objective
        Xm = X0; Ym = Y0;
        Yt = Y0'; Ys = Y0*Yt; MYt = M*Yt;
    else
        % extrapolation
        w = (t0-1)/t; % extrapolation weight
        wx = min([w,rw*sqrt(Lx0/Lx)]); % choose smaller one for convergence
        wy = min([w,rw*sqrt(Ly0/Ly)]);
        Xm = X + wx*(X-X0);   Ym = Y + wy*(Y-Y0); % extrapolation
        X0 = X; Y0 = Y; t0 = t; obj0 = obj;
    end 
end
out.iter = k;
fprintf('\n');  % report # of iterations

