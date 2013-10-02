function h = h_estimate_Gustafsson_exp(Yerr, h_old, ewt, p, hflag)
% Usage: h = h_estimate_Gustafsson_exp(Yerr, h_old, ewt, p, hflag)
%
% Adaptive time step estimation routine, that attempts to guarantee a 
% local truncation error satisfying the bound
%    norm(local error,inf) < norm(rtol*y_i + atol_i,inf)
%
% Inputs:
%        Yerr -- estimated error in time-evolved solution (vector)
%       h_old -- previous time step size (scalar)
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
   clear r_n;
   h = 1;
   return;
end

% set variables
dt_growth = 10;
kS = 1.0/p;
kI = 0.3/p;
kP = 0.4/p;
persistent r_n;

% first step
if isempty(r_n)

   % compute new error estimate, time step
   r_n = max(norm(Yerr.*ewt,inf), eps);
   h = r_n^(-kS) * h_old;

% general step
else

   % store error estimate from previous step
   r_save = r_m;
   r_m = r_n;
   
   % compute new error estimate
   r_n = max(norm(Yerr.*ewt,inf), eps);

   % compute new time step
   h = r_n^(-kI) * (r_n/r_m)^(-kP) * h_old;

   % revert to original history if requested
   if (hflag == -1)
      r_n = r_m;
      r_m = r_save;
   end
   
end

% enforce growth factor
h = min(dt_growth*h_old, h);


% end of function