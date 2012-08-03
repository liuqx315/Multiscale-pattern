function [tvals,Y,nsteps] = solve_DIRK(fcn, Jfcn, tvals, Y0, B, rtol, atol, hmin, hmax, hmethod)
% usage: [tvals,Y,nsteps] = solve_DIRK(fcn, Jfcn, tvals, Y0, B, rtol, atol, hmin, hmax, hmethod)
%
% Adaptive time step DIRK solver for the vector-valued ODE problem
%     y' = F(t,Y), t in tspan,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:  fcn = function name for ODE right-hand side, F(t,Y)
%          Jfcn = function name for Jacobian of ODE right-hand side, J(t,Y)
%          tvals = [t0, t1, t2, ..., tN]
%          Y0 = initial values
%          B = Butcher matrix for IRK coefficients, of the form
%                 B = [c A;
%                      0 b;
%                      0 b2 ]
%              The b2 row is optional, and provides coefficients for an
%              embedded method.
%          rtol = desired time accuracy relative tolerance 
%          atol = desired time accuracy absolute tolerance 
%          hmin = min internal time step size (must be smaller than t(i)-t(i-1))
%          hmax = max internal time step size (can be smaller than t(i)-t(i-1))
%          hmethod = integer flag denoting which time adaptivity strategy to use
%
% Outputs: t = tspan
%          y = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%              y(t*) is a column vector of length m.
%          nsteps = number of internal time steps taken
%
% Note: to run in fixed-step mode, just call the solver using hmin=hmax as
% the desired time step size.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% get number of stages for DIRK method
[Brows, Bcols] = size(B);
s = Bcols - 1;
if (Brows > Bcols) 
   embedded = 1;
else
   embedded = 0;
end

% extract DIRK method information from B
c = B(1:s,1);
b = (B(s+1,2:s+1))';
A = B(1:s,2:s+1);

% estimate method order (and embedding order)
if (s == 4) 
   q_method = 3;
   p_method = 2;
elseif (s == 5)
   q_method = 4;
   p_method = 3;
elseif (s == 6)
   q_method = 4;
   p_method = 3;
elseif (s == 7)
   q_method = 5;
   p_method = 4;
elseif (s == 8)
   q_method = 5;
   p_method = 4;
else
   q_method = 5;
   p_method = 4;
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
newt_maxit = 20;
%newt_tol   = 1e-10;
%newt_tol   = 1e-12;
newt_tol   = 1e-11;
newt_alpha = 1;
dt_reduce  = 0.1;
ONEMSM     = 1.0-sqrt(eps);
ONEPSM     = 1.0+sqrt(eps);
ERRTOL     = 1.2;

% store temporary states
t = tvals(1);
Ynew = Y0;

% create Fdata structure
Fdata.fname = fcn;
Fdata.Jname = Jfcn;
Fdata.B = B;
Fdata.s = s;

% set initial time step size
h = hmin;

% reset time step controller
h_estimate(0, 0, 0, 0, 0, 0, hmethod, 1);


% initialize work counter
nsteps = 0;

% iterate over time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep)*ONEMSM)
      
      % limit internal time step if needed
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time
      Fdata.h = h;
      Fdata.yold = Y0;

      % initialize data storage for multiple stages
      z = zeros(m,s);

      % reset stage failure flag
      st_fail = 0;

% $$$       fprintf('\n');
% $$$       fprintf('  Attempting internal step of size h = %g\n',h);      
      
      % loop over stages
      for stage=1:s
	 
	 % set stage initial guess as previous stage solution
	 Yguess = Ynew;
      
	 % set stage number into Fdata structure
	 Fdata.stage = stage;
	 
	 % construct RHS comprised of old time data
	 %    zi = y_n + h*sum_{j=1}^s (a(i,j)*fj)
	 %    zi = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj) + h*(a(i,i)*fi)
	 %    zi - h*(a(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
	 % =>
	 %    rhs = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
	 Fdata.rhs = Y0;
	 for j=1:stage-1
	    Fdata.rhs = Fdata.rhs + h*A(stage,j)*feval(fcn, t+h*c(j), z(:,j));
	 end

% $$$          fprintf('    stage %i: yguess = ', stage);
% $$$          for entry=1:length(Yguess), fprintf('%g, ',Yguess(entry)); end
% $$$          fprintf('\n');

% $$$          fprintf('             rhs = ');
% $$$          for entry=1:length(Fdata.rhs), fprintf('%g, ',Fdata.rhs(entry)); end
% $$$          fprintf('\n');


	 % call newton solver to compute new stage solution
	 Fdata.t = t;
	 [Ynew,inewt,ierr] = newton_damped('F_DIRK', 'A_DIRK', Yguess, ...
	     Fdata, newt_tol, newt_maxit, newt_alpha);

% $$$          fprintf('             ynew = ');
% $$$          for entry=1:length(Ynew), fprintf('%g, ',Ynew(entry)); end
% $$$          fprintf('\n');

	 % check newton error flag, if failure break out of stage loop
	 if (ierr ~= 0) 
	    % fprintf('solve_DIRK warning: stage failure, cutting timestep\n');
	    st_fail = 1;
	    c_fails = c_fails + 1;
	    break;
	 end
	 
	 % update stored solution with new value
	 z(:,stage) = Ynew;
	 
      end
      nsteps = nsteps + 1;
      
      
      % compute new solution, embedding (if available)
      [Ynew,Y2] = Y_DIRK(z,Fdata);

% $$$       fprintf('  step solution: Ynew = ');
% $$$       for entry=1:length(Ynew), fprintf('%g, ',Ynew(entry)); end
% $$$       fprintf('\n');

      % check whether step accuracy meets tolerances (only if stages successful)
      if ((st_fail == 0) & embedded)

	 % compute error in current step
	 err_step = max(norm((Ynew - Y2)./(rtol*Ynew + atol),inf), eps);
	 
% $$$          fprintf('  error estimate = %g\n',err_step);

	 % if error too high, flag step as a failure (to be recomputed)
	 if (err_step > ERRTOL*ONEPSM) 
	    a_fails = a_fails + 1;
	    st_fail = 1;
	 end
	 
      end


      % if step was successful
      if (st_fail == 0) 
      
	 % update old solution, current time
	 Y0 = Ynew;
	 t = t + h;
   
	 % estimate error and update time step
	 if (embedded) 
	    h = h_estimate(Ynew, Y2, h, rtol, atol, q_method, hmethod, 0);
	 else
	    h = hmin;
	 end

% $$$          fprintf('  successful step, next h = %g\n',h);

      % if step failed, reduce step size and retry
      else
	 
	 % reset solution guess, update work counter, reduce time step
	 Ynew = Y0;
	 h = h * dt_reduce;
         if (h <= hmin*ONEMSM) 
            return
         end

% $$$          fprintf('    failed step, next h = %g\n',h);

      end
      
   end

% $$$    fprintf('Output solution = ');
% $$$    for entry=1:length(Ynew), fprintf('%g, ',Ynew(entry)); end
% $$$    fprintf('\n');

   % store updated solution
   Y(:,tstep) = Ynew;
   
end

% $$$ fprintf('solve_DIRK: step failures:  %i convergence,  %i accuracy\n',...
% $$$     c_fails,a_fails);

% end function
