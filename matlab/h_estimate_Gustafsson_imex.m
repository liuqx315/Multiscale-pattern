function h = h_estimate_Gustafsson_imex(Yerr, h_old, ewt, p, hflag)
% Usage: h = h_estimate_Gustafsson_imex(Yerr, h_old, ewt, p, hflag)
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


% get estimates from both implicit and explicit approaches
h_imp = h_estimate_Gustafsson(Yerr, h_old, ewt, p, hflag);
h_exp = h_estimate_Gustafsson_exp(Yerr, h_old, ewt, p, hflag);

% return minimum of these steps
h = min([h_imp, h_exp]);

% end of function