function [tvals,Y,nsteps] = solve_DIRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,hmethod)
% usage: [tvals,Y,nsteps] = solve_DIRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,hmethod)
%
% Adaptive time step diagonally-implicit Runge-Kutta solver for the
% vector-valued ODE problem 
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
%     nsteps = number of internal time steps taken by method
%
% Note: to run in fixed-step mode, call with hmin=hmax as the desired 
% time step size, and set the tolerances to large positive numbers.
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

DO_OUTPUT = 0;

% extract DIRK method information from B
[Brows, Bcols] = size(B);
s = Bcols - 1;        % number of stages
c = B(1:s,1);         % stage time fraction array
b = (B(s+1,2:s+1))';  % solution weights (convert to column)
A = B(1:s,2:s+1);     % RK coefficients

% initialize as non-embedded, until proven otherwise
embedded = 0;
p = 0;
if (Brows > Bcols)
   if (max(abs(B(s+2,2:s+1))) > eps)
      embedded = 1;
      b2 = (B(s+2,2:s+1))';
      p = B(s+2,1);
   end
end

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% initialize diagnostics
c_fails = 0;
a_fails = 0;

% set the solver parameters
newt_maxit = 20;             % max number of Newton iterations
newt_tol   = 1e-10;          % Newton solver tolerance
newt_alpha = 1;              % Newton damping parameter
h_reduce   = 0.1;            % failed step reduction factor 
SMALL      = sqrt(eps);      % tolerance for floating-point comparisons
ONEMSM     = 1.0-SMALL;      % coefficients to account for
ONEPSM     = 1.0+SMALL;      %   floating-point roundoff
ERRTOL     = 1.1;            % upper bound on allowed step error
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
ewt = 1.0/(rtol*Ynew + atol);

if (DO_OUTPUT)
   fprintf('updated ewt =');
   for entry=1:length(ewt), fprintf('%19.16g, ',ewt(entry)); end
   fprintf('\n');
end

% reset time step controller
h_estimate(0, 0, 0, 0, hmethod, 1);

% initialize work counter
nsteps = 0;

% initialize error failure counters
small_nef = 2;
nef = 0;
last_fail = 0;

% iterate over time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep)*ONEMSM)
      
      % bound internal time step 
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time

      % look-ahead to avoid a small step to output time
      if (abs(tvals(tstep)-(t+h)/tvals(tstep)) < SMALL)
         h = tvals(tstep)-t;
      end
      
      % set Fdata values for this step
      Fdata.h    = h;    % current step size
      Fdata.yold = Y0;   % solution from previous step
      Fdata.t    = t;    % time of last successful step

      % initialize data storage for multiple stages
      z = zeros(m,s);

      % reset stage failure flag
      st_fail = 0;
      
      if (DO_OUTPUT)
         fprintf('\n');
         fprintf('  Attempting internal step, t = %19.16g, h = %19.16g\n',t,h);
      end
      
      % loop over stages
      for stage=1:s
	 
         % set Newton initial guess as previous stage solution
         Yguess = Ynew;
         
         % set current stage index into Fdata structure
         Fdata.stage = stage;
         
         % construct RHS comprised of old time data
         %    zi = y_n + h*sum_{j=1}^s (a(i,j)*fj)
         % <=>
         %    zi - h*(a(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
         % =>
         %    rhs = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
         Fdata.rhs = Y0;
         for j = 1:stage-1
            Fdata.rhs = Fdata.rhs + h*A(stage,j)*feval(fcn, t+h*c(j), z(:,j));
         end
         
         if (DO_OUTPUT)
            fprintf('    stage %i: tstage = %19.16g\n', stage, t+h*c(stage));
            fprintf('         yguess = ');
            for entry=1:length(Yguess), fprintf('%19.16g, ',Yguess(entry)); end
            fprintf('\n');
         end

         % call Newton solver to compute new stage solution
         [Ynew,lin,ierr] = newton_damped('F_DIRK', 'A_DIRK', Yguess, Fdata, ...
                                         newt_tol, newt_maxit, newt_alpha);

         if (DO_OUTPUT)
            fprintf('           ynew = ');
            for entry=1:length(Ynew), fprintf('%19.16g, ',Ynew(entry)); end
            fprintf('\n');
         end
         
         % if Newton method failed, set relevant flags/statistics
         % and break out of stage loop
         if (ierr ~= 0) 
            st_fail = 1;
            c_fails = c_fails + 1;
            break;
         end
         
         % store stage solution
         z(:,stage) = Ynew;
         
      end  % end stage loop

      % increment number of internal time steps taken
      nsteps = nsteps + 1;
       
      % compute new solution (and embedding if available)
      [Ynew,Yerr] = Y_DIRK(z,Fdata);

      if (DO_OUTPUT)
         fprintf('  step solution: Ynew = ');
         for entry=1:length(Ynew), fprintf('%19.16g, ',Ynew(entry)); end
         fprintf('\n');
      end
 
      % if stages succeeded and time step adaptivity enabled, check step accuracy
      if ((st_fail == 0) & embedded)

	 % compute error in current step
	 err_step = norm(Yerr.*ewt,inf);
	 
         if (DO_OUTPUT)
            fprintf('  error estimate = %19.16g\n',err_step);
         end
         
         % if error too high, flag step as a failure (will be be recomputed)
         if (err_step > ERRTOL*ONEPSM) 
            a_fails = a_fails + 1;
            st_fail = 1;
         end
         
      end

      % if step was successful (solves succeeded, and error acceptable)
      if (st_fail == 0) 

         % reset error failure counter
         nef = 0;
         
         % update solution and time for last successful step
         Y0 = Ynew;
         t  = t + h;

	 % estimate error and update time step
	 if (embedded) 
	    hnew = h_estimate(Yerr, h, ewt, p, hmethod, 0);
	 else
	    hnew = hmin;
         end
  
         % if last step failed, disallow growth on this step and
         % turn off flag, otherwise just update as normal
         if (last_fail)
            h = min(hnew, h);
            last_fail = 0;
         else
            h = hnew;
         end

         % update error weight vector
         ewt = 1.0/(rtol*Ynew + atol);
         
         if (DO_OUTPUT)
            fprintf('updated ewt =');
            for entry=1:length(ewt), fprintf('%19.16g, ',ewt(entry)); end
            fprintf('\n');

            fprintf('  successful step, next h = %19.16g\n',h);
         end

      % if step solves or error test failed
      else
	 
         % if already at minimum step, just return with failure
         if (h <= hmin) 
            fprintf('Cannot achieve desired accuracy.\n');
            fprintf('Consider reducing hmin or increasing rtol.\n');
            return
         end

         % update error failure counter
         nef = nef + 1;
         last_fail = 1;
         
         % update time step
         if (nef >= small_nef)
            h = h * h_reduce;
         else
	    h = h_estimate(Yerr, h, ewt, p, hmethod, -1);
         end

         % reset guess, and try solve again
         Ynew = Y0;

      end  % end logic tests for step success/failure
      
   end  % end while loop attempting to solve steps to next output time

   if (DO_OUTPUT)
      fprintf('Output solution = ');
      for entry=1:length(Ynew), fprintf('%19.16g, ',Ynew(entry)); end
      fprintf('\n');
   end

   % store updated solution in output array
   Y(:,tstep) = Ynew;
   
end  % time step loop

% $$$ fprintf('solve_DIRK: step failures:  %i convergence,  %i accuracy\n',...
% $$$     c_fails,a_fails);

% end function
