function h = h_estimate_Gustafsson_exp(Y1, Y2, h_old, r_tol, a_tol, p)
% Usage: h = h_estimate_Gustafsson_exp(Y1, Y2, h_old, r_tol, a_tol, p)
%
% Adaptive time step estimation routine, that attempts to guarantee a 
% local truncation error satisfying the bound
%    norm(local error,inf) < norm(rtol*y_i + atol_i,inf)
%
% Inputs:  Y1 -- first time-evolved solution (more accurate)
%          Y2 -- second time-evolved solution (same size as Y1)
%       h_old -- previous time step size (scalar)
%       r_tol -- desired relative tolerance (scalar)
%       a_tol -- desired absolute tolerance (vector, size of Y1)
%           p -- order of accuracy for predictor
%
% Output:   h -- new time step
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2011
% All Rights Reserved


% set variables
dt_growth = 10;
kS = 1.0/p;
kI = 0.3/p;
kP = 0.4/p;
persistent r_n;

% first step
if isempty(r_n)

   % compute new error estimate, time step
   r_n = max(norm((Y1 - Y2)./(r_tol*Y1 + a_tol),inf), eps);
   h = r_n^(-kS) * h_old;

% general step
else

   % store error estimate from previous step
   r_m = r_n;
   
   % compute new error estimate
   r_n = max(norm((Y1 - Y2)./(r_tol*Y1 + a_tol),inf), eps);

   % compute new time step
   h = r_n^(-kI) * (r_n/r_m)^(-kP) * h_old;
   
end

% enforce growth factor
h = min(dt_growth*h_old, h);


% end of function