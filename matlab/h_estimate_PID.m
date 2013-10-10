function h = h_estimate_PID(Yerr, h_old, ewt, p, hflag)
% Usage: h = h_estimate_PID(Yerr, h_old, ewt, p, hflag)
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
   clear Eratio Eratio_old;
   h = 1;
   return;
end

% set variables
safety = 0.9;
dt_growth = 10;
alpha = 0.49/p;
beta  = 0.34/p;
gamma = 0.1/p;
persistent Eratio;
if isempty(Eratio)
   Eratio = 1;
end
persistent Eratio_old;
if isempty(Eratio_old)
   Eratio_old = 1;
end

% update old error estimates
Eratio_oldest = Eratio_older;
Eratio_older = Eratio_old;
Eratio_old = Eratio;

% compute new error estimate ratio, bound from below 
Eratio = max(norm(Yerr.*ewt,inf), eps);

% compute updated time step
h = safety * h_old * Eratio^(-alpha) * Eratio_old^(beta) * Eratio_older^(-gamma);

% enforce growth factor
h = min(dt_growth*h_old, h);

% revert to original history if requested
if (hflag == -1)
   Eratio = Eratio_old;
   Eratio_old = Eratio_older;
   Eratio_older = Eratio_oldest;
end


% end of function