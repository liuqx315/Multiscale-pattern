function h = h_estimate_Gustafsson(Y1, Y2, h_cur, r_tol, a_tol, p, reset)
% Usage: h = h_estimate_Gustafsson(Y1, Y2, h_cur, r_tol, a_tol, p, reset)
%
% Adaptive time step estimation routine, that attempts to guarantee a 
% local truncation error satisfying the bound
%    norm(local error,inf) < norm(rtol*y_i + atol_i,inf)
%
% Inputs:  Y1 -- first time-evolved solution (more accurate)
%          Y2 -- second time-evolved solution (same size as Y1)
%       h_cur -- previous time step size (scalar)
%       r_tol -- desired relative tolerance (scalar)
%       a_tol -- desired absolute tolerance (vector, size of Y1)
%           p -- order of accuracy for predictor
%       reset -- flag to denote reset of history
%
% Output:   h -- new time step
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2011
% All Rights Reserved


if (reset == 1)
   clear r_p;
   clear h_old;
   h = 1;
   return;
end

% set variables
dt_growth = 10;
k1 = 1.0/p;
k2 = 1.0/p;
persistent r_p;
persistent h_old;

% first step
if isempty(r_p)

   % compute new error estimate, time step
   r_p = max(norm((Y1 - Y2)./(r_tol*Y1 + a_tol),inf), eps);
   h = r_p^(-k2) * h_cur;

% general step
else

   % store error estimate from previous step
   r_n = r_p;
   
   % compute new error estimate
   r_p = max(norm((Y1 - Y2)./(r_tol*Y1 + a_tol),inf), eps);

   % compute new time step
   h = (h_cur/h_old) * r_p^(-k2) * (r_p/r_n)^(-k1) * h_cur;
   
end

% store current time step for later
h_old = h_cur;

% enforce growth factor
h = min(dt_growth*h_old, h);


% end of function