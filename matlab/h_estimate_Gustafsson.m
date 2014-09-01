function h = h_estimate_Gustafsson(Yerr, h_cur, ewt, p, hflag)
% Usage: h = h_estimate_Gustafsson(Yerr, h_cur, ewt, p, hflag)
%
% Adaptive time step estimation routine, that attempts to guarantee a 
% local truncation error satisfying the bound
%    norm(local error,inf) < norm(rtol*y_i + atol_i,inf)
%
% Inputs:
%        Yerr -- estimated error in time-evolved solution (vector)
%       h_cur -- previous time step size (scalar)
%         ewt -- error weight vector (encodes tolerances)
%           p -- order of accuracy for predictor
%       hflag -- flag to denote what to do with history:
%                  1 => reset,  -1 => leave alone,  0 => update
%
% Output:   h -- new time step
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------


if (hflag == 1)
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
   r_p = max(norm(Yerr.*ewt,inf), eps);
   h = r_p^(-k2) * h_cur;

% general step
else

   % store error estimate from previous step
   r_n = r_p;
   
   % compute new error estimate
   r_p = max(norm(Yerr.*ewt,inf), eps);

   % compute new time step
   h = (h_cur/h_old) * r_p^(-k2) * (r_p/r_n)^(-k1) * h_cur;
   
   % revert to original history if requested
   if (hflag == -1)
      r_p = r_n;
   end
   
end

% store current time step for later
h_save = h_old;
h_old = h_cur;

% enforce growth factor
h = min(dt_growth*h_old, h);

% revert to original history if requested
if (hflag == -1)
   h_old = h_save;
end
   

% end of function