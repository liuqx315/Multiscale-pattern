function [tvals,Y,nsteps] = solve_IRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,hmethod)
% usage: [tvals,Y,nsteps] = solve_IRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,hmethod)
%
% Fixed time step implicit Runge-Kutta solver for the vector-valued
% ODE problem 
%     y' = F(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fcn    = string holding function name for F(t,Y)
%     Jfcn   = string holding function name for Jacobian of F, J(t,Y)
%     tvals  = [t0, t1, t2, ..., tN]
%     Y0     = initial value array (column vector of length m)
%     B      = Butcher matrix for IRK coefficients, of the form
%                 B = [c A;
%                      q b;
%                      p b2 ]
%              Here, c is a vector of stage time fractions (s-by-1),
%                    A is a matrix of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    b is a vector of solution weights (1-by-s),
%                    p is an integer denoting the embedding order of accuracy,
%                    b2 is a vector of embedding weights (1-by-s),
%              The [p, b2] row is optional.  If that row is not
%              provided the method will default to taking fixed
%              step sizes of size hmin.
%     rtol   = desired relative error of solution  (scalar)
%     atol   = desired absolute error of solution  (vector or scalar)
%     hmin   = minimum internal time step size (hmin <= t(i)-t(i-1), for all i)
%     hmax   = maximum internal time step size (hmax >= hmin)
%    hmethod = integer flag denoting which time adaptivity strategy to use
%
% Outputs: 
%     tvals  = the same as the input array tvals
%     y      = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%               y(t*) is a column vector of length m.
%     nsteps = number of internal time steps taken by method.
%              Note: if adaptivity is enabled, this includes the
%              two half steps used in the Richardson error estimate
%              and extrapolation.
%
% Note: since embeddings in IRK methods are typically *much* lower
% order than the original method, we instead base our adaptivity on
% Richardson extrapolation.  Specifically, for every time step, we
% solve using both a single step of size h and two steps of size
% h/2.  The difference of these provides an estimate on the local
% error,  Moreover, an appropriately-chosen linear combination of
% these provides a solution that is one order of accuracy higher
% than the method itself -- this is the solution that is stored and
% returned to the user.
%
% Note2: to run in fixed-step mode, call with hmin=hmax as the desired 
% time step size, and set the tolerances to large positive numbers.
%   
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

   
% get number of stages and internal time fractions for IRK method
[Brows, Bcols] = size(B);
s = Bcols - 1;        
c = B(1:s,1);

% check whether time step adaptivity is desired
adaptivity = 0;
if (abs(hmax-hmin)/abs(hmax) > sqrt(eps))
   q = B(s+1,1);     % order of accuracy for method
   c1 = -1/(2^q-1);  % Richardson extrapolation factor for h step
   c2 = 1 - c1;      % Richardson extrapolation factor for h/2 step
   adaptivity = 1;
end

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% initialize diagnostics
c_fails = 0;   % total convergence failures
a_fails = 0;   % total accuracy failures

% set the solver parameters
newt_maxit = 20;           % max number of Newton iterations
newt_tol   = 1e-10;        % Newton solver tolerance
newt_alpha = 1;            % Newton damping parameter
h_reduce   = 0.1;          % failed step reduction factor 
ONEMSM     = 1-sqrt(eps);  % coefficients to account for
ONEPSM     = 1+sqrt(eps);  %   floating-point roundoff
ERRTOL     = 1.1;          % upper bound on allowed step error
                           %   (in WRMS norm)

% initialize temporary variables
t = tvals(1);
Ynew = Y0;

% create Fdata structure for Newton solver and step solutions
Fdata.fname = fcn;    % ODE RHS function name
Fdata.Jname = Jfcn;   % ODE RHS Jacobian function name
Fdata.B     = B;      % Butcher table 
Fdata.s     = s;      % number of stages

% set initial time step size
h = hmin;

% initialize error weight vector
ewt = 1.0./(rtol*Ynew + atol);

% reset time step controller
h_estimate(0, 0, 0, 0, hmethod, 1);

% initialize work counter
nsteps = 0;

% iterate over time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep)*ONEMSM)
      
      % bound internal time step 
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time

      % set Fdata values for this step
      Fdata.h    = h;      % current step size
      Fdata.yold = Ynew;   % solution from previous step
      Fdata.t    = t;      % time of last successful step

      % reset solve failure flag
      st_fail = 0;
      
      % solve with time step h
      
      % set Newton initial guesses as previous step solution
      z = zeros(s*m,1);
      for i = 0:s-1
	 z(i*m+1:(i+1)*m) = Ynew;
      end

      % call Newton solver to update solution in time
      [z,lin,ierr] = newton_damped('F_IRK', 'A_IRK', z, Fdata, ...
                                   newt_tol, newt_maxit, newt_alpha);

      % increment step statistics 
      nsteps = nsteps + 1;
      
      % compute solution with this h
      Y1 = Y_IRK(z,Fdata);

      % if Newton method failed, set relevant flags/statistics
      if (ierr ~= 0) 
         st_fail = 1;
         c_fails = c_fails + 1;
      end
         
      % if adaptivity enabled, solve with steps of size h/2
      if (adaptivity & (st_fail == 0)) 
         
         % set Fdata values for this half-step
         Fdata.h    = 0.5*h;  % half-step size
         Fdata.yold = Ynew;   % solution from previous step
         Fdata.t    = t;      % time of last successful step
      
         % set Newton initial guesses as linear interpolants of
         % full step solutions
         z = zeros(s*m,1);
         for i = 0:s-1
            ti = c(i+1)*0.5;
            z(i*m+1:(i+1)*m) = (1-ti)*Ynew + ti*Y1;
         end
         
         % call Newton solver to update solution in time
         [z,lin,ierr] = newton_damped('F_IRK', 'A_IRK', z, Fdata, ...
                                      newt_tol, newt_maxit, newt_alpha);
         
         % increment step statistics
         nsteps = nsteps + 1;
         
         % compute half-step solution
         Y2 = Y_IRK(z,Fdata);
         
         % if Newton method failed, set relevant flags/statistics
         if (ierr ~= 0) 
            st_fail = 1;
            c_fails = c_fails + 1;
         end

         % if first half-step succeeded, take second half-step
         if (st_fail == 0) 
            
            % set Fdata values for second half-step
            Fdata.h    = 0.5*h;    % half-step size
            Fdata.yold = Y2;       % solution from previous half-step
            Fdata.t    = t+0.5*h;  % time of half-step 
            
            % set Newton initial guesses as linear interpolants of
            % full step solutions
            z = zeros(s*m,1);
            for i = 0:s-1
               ti = 0.5 + c(i+1)*0.5;
               z(i*m+1:(i+1)*m) = (1-ti)*Ynew + ti*Y1;
            end
            
            % call Newton solver to update solution in time
            [z,lin,ierr] = newton_damped('F_IRK', 'A_IRK', z, Fdata, ...
                                         newt_tol, newt_maxit, newt_alpha);
         
            % increment step statistics
            nsteps = nsteps + 1;
         
            % compute full-step solution (store back in Y2)
            Y2 = Y_IRK(z,Fdata);
         
            % if Newton method failed, set relevant flags/statistics
            if (ierr ~= 0) 
               st_fail = 1;
               c_fails = c_fails + 1;
            end

         end % end second half-step
         
      end % half-step solutions
      
      % if solves succeeded and time step adaptivity enabled, check step accuracy
      if (adaptivity & (st_fail == 0))

         % estimate error in current step
         Yerr = Y1-Y2;
         err_step = max(norm(Yerr.*ewt,inf), eps);
         
         % if error too high, flag step as a failure (will be be recomputed)
         if (err_step > ERRTOL*ONEPSM) 
            a_fails = a_fails + 1;
            st_fail = 1;
         end
         
      end

      % if step was successful (solves succeeded, and error acceptable)
      if (st_fail == 0) 
         
         % update time for last successful step
         t  = t + h;
         
         % update solution (use Richardson extrapolation if available)
         if (adaptivity)
            Ynew = c1*Y1 + c2*Y2;
         else
            Ynew = Y1;
         end
         Y0 = Ynew;
         
         % update error weight vector
         ewt = 1.0./(rtol*Ynew + atol);
         
	 % estimate error and update time step
         if (adaptivity) 
	    h = h_estimate(Yerr, h, ewt, p, hmethod, 0);
         % otherwise, just use the fixed minimum input step size
         else
            h = hmin;
         end
         
      % if step solves or error test failed
      else

         % if already at minimum step, just return with failure
         if (h <= hmin) 
            fprintf('Cannot achieve desired accuracy.\n');
            fprintf('Consider reducing hmin or increasing rtol.\n');
            return
         end

         % otherwise, reset guess, reduce time step, retry solve
         Ynew = Y0;
         h = h * h_reduce;
         
      end  % end logic tests for step success/failure
      
   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep) = Ynew;
   
end  % time step loop


% end solve_IRK function
end





function y = Y_IRK(z, Fdata)
% usage: [y,y2] = Y_IRK(z, Fdata)
%
% Inputs:
%    z     = stage solutions [z1, ..., zs]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    y     = step solution built from the z values
%    y2    = embedded solution (if embedding included in Butcher 
%               table; otherwise the same as y)

% extract method information from Fdata
B = Fdata.B;
[Brows, Bcols] = size(B);
s = Bcols - 1;
c = B(1:s,1);
b = (B(s+1,2:s+1))';
A = B(1:s,2:s+1);

% get some problem information
zlen = length(z);
nvar = floor(zlen/s);
if (nvar*s ~= zlen)
   error('Y_IRK error: input has incorrect length (must be a multiple of s)');
end

% reshape our z arguments into separate vectors for each stage
z = reshape(z,nvar,s);

% call f at our stages
f = zeros(nvar,s);
for is=1:s
   t = Fdata.t + Fdata.h*c(is);
   f(:,is) = feval(Fdata.fname, t, z(:,is));
end

% form the solution
%    ynew = yold + h*sum(b(j)*fj)
y = Fdata.yold + Fdata.h*f*b;

% end of function
end

