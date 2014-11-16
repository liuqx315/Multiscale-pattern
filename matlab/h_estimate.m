function h = h_estimate(Yerr, h_old, ewt, p, hmethod, hflag)
% Usage: h = h_estimate(Yerr, h_old, ewt, p, hmethod, hflag)
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
%     hmethod -- adaptivity strategy to use
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

% call relevant method
switch (hmethod)
   
  case(1)  % I controller
    h = h_estimate_I(Yerr, h_old, ewt, p, hflag);
    
  case(2)  % PI controller
    h = h_estimate_PI(Yerr, h_old, ewt, p, hflag);
    
  case(3)  % PID controller
    h = h_estimate_PID(Yerr, h_old, ewt, p, hflag);
    
  case(4)  % implicit Gustafsson controller
    h = h_estimate_Gustafsson(Yerr, h_old, ewt, p, hflag);
    
  case(5)  % explicit Gustafsson controller
    h = h_estimate_Gustafsson_exp(Yerr, h_old, ewt, p, hflag);

  case(6)  % imex Gustafsson controller
    h = h_estimate_Gustafsson_imex(Yerr, h_old, ewt, p, hflag);
    
end
    

% end of function