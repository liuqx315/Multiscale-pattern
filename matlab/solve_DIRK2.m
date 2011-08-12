function [tvals,Y,nsteps] = solve_DIRK2(fcn, Jfcn, tvals, Y0, B, tol, hmin, hmax)
% usage: [tvals,Y,nsteps] = solve_DIRK2(fcn, Jfcn, tvals, Y0, B, tol, hmin, hmax)
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
%          tol = desired time accuracy tolerance 
%          hmin = min internal time step size (must be smaller than t(i)-t(i-1))
%          hmax = max internal time step size (can be smaller than t(i)-t(i-1))
%          nsteps = number of internal time steps taken (total stage steps)
%
% Outputs: t = tspan
%          y = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%              y(t*) is a column vector of length m.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% get number of stages for IRK method
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

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% set the solver parameters
newt_maxit = 20;
newt_tol   = 1e-10;
newt_alpha = 1;
dt_safety  = 0.1;
dt_reduce  = 0.1;

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

% initialize work counter
nsteps = 0;

% iterate over time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep))
      
      % limit internal time step if needed
      h = min([h, hmax, tvals(tstep)-t]);
      Fdata.h = h;
      Fdata.yold = Y0;

      % initialize data storage for multiple stages
      z = zeros(m,s);

      % reset stage success flag
      st_success = 0;
      
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
	 
	 % call newton solver to compute new stage solution
	 Fdata.t = t;
	 [Ynew,ierr] = newton_damped('F_DIRK', 'A_DIRK', Yguess, Fdata, ...
	     newt_tol, newt_maxit, newt_alpha);
	 nsteps = nsteps + 1;
	 
	 % check newton error flag, if failure break out of stage loop
	 if (ierr ~= 0) 
	    fprintf('solve_DIRK2 warning: stage failure, cutting timestep\n');
	    st_success = 1;
	    break;
	 end
	 
	 % update stored solution with new value
	 z(:,stage) = Ynew;
	 
      end

      % check whether all stages solved successfully
      if (st_success == 0) 
      
	 % use stage solutions to determine time-evolved answer, update old
	 [Ynew,Y2] = Y_DIRK(z,Fdata);
	 Y0 = Ynew;

	 % update time, work counter
	 t = t + h;
   
	 % estimate error and update time step, assuming that method order
	 % equals number of stages-1 (this should be an input argument)
	 if (embedded) 
	    err = norm(Y2-Ynew,inf) + sqrt(eps)*tol;
	    h = max([h * (dt_safety*tol/err)^(1/(s-1)), hmin]);
	 end
	 
      % if any stages failed, reduce step size and retry
      else
	 
	 % reset solution guess, update work counter, reduce time step
	 Ynew = Y0;
	 h = h*dt_reduce;
      
      end
      
   end

   % store updated solution
   Y(:,tstep) = Ynew;
   
end


% end function
