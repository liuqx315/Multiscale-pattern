function h = h_estimate(Y1, Y2, h_old, r_tol, a_tol, p, hmethod, reset)
% Usage: h = h_estimate(Y1, Y2, h_old, r_tol, a_tol, p, hmethod, reset)
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
%     hmethod -- adaptivity strategy to use
%       reset -- flag to denote reset of history
%
% Output:   h -- new time step
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2011
% All Rights Reserved

% call relevant method
switch (hmethod)
   
  case(1)  % I controller
    h = h_estimate_I(Y1, Y2, h_old, r_tol, a_tol, p, reset);
    
  case(2)  % PI controller
    h = h_estimate_PI(Y1, Y2, h_old, r_tol, a_tol, p, reset);
    
  case(3)  % PID controller
    h = h_estimate_PID(Y1, Y2, h_old, r_tol, a_tol, p, reset);
    
  case(4)  % implicit Gustafsson controller
    h = h_estimate_Gustafsson(Y1, Y2, h_old, r_tol, a_tol, p, reset);
    
  case(5)  % explicit Gustafsson controller
    h = h_estimate_Gustafsson_exp(Y1, Y2, h_old, r_tol, a_tol, p, reset);

  case(6)  % imex Gustafsson controller
    h = h_estimate_Gustafsson_imex(Y1, Y2, h_old, r_tol, a_tol, p, reset);
    
end
    

% end of function