function h = h_estimate_Gustafsson_imex(Y1, Y2, h_old, r_tol, a_tol, p, reset)
% Usage: h = h_estimate_Gustafsson_imex(Y1, Y2, h_old, r_tol, a_tol, p, reset)
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
%       reset -- flag to denote reset of history
%
% Output:   h -- new time step
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2011
% All Rights Reserved


% get estimates from both implicit and explicit approaches
h_imp = h_estimate_Gustafsson(Y1, Y2, h_old, r_tol, a_tol, p, reset);
h_exp = h_estimate_Gustafsson_exp(Y1, Y2, h_old, r_tol, a_tol, p, reset);

% return minimum of these steps
h = min([h_imp, h_exp]);

% end of function