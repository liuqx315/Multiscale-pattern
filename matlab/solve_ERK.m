function [tvals,Y,nsteps] = solve_ERK(fcn, EStabFn, tvals, Y0, B, rtol, atol, hmin, hmax)
% usage: [tvals,Y,nsteps] = solve_ERK(fcn, EStabFn, tvals, Y0, B, rtol, atol, hmin, hmax)
%
% Adaptive time step explicit RK solver for the vector-valued ODE problem
%     y' = F(t,Y), t in tspan,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:  fcn = function name for ODE right-hand side,
%          EStabFn = function name for stability constraint on fcn
%          tvals = [t0, t1, t2, ..., tN]
%          Y0 = initial values
%          B = Butcher matrix for ERK coefficients, of the form
%                 B = [c A;
%                      0 b;
%                      0 b2 ]
%              The b2 row is optional, but if not provided the method will
%              default to taking fixed step sizes, at min(hmin, hstab),
%              where hstab is the result of EStabFn(t,Y). 
%          rtol = desired time accuracy relative tolerance 
%          atol = desired time accuracy absolute tolerance 
%          hmin = min internal time step size (must be smaller than t(i)-t(i-1))
%          hmax = max internal time step size (can be smaller than t(i)-t(i-1))
%          nsteps = number of internal time steps taken (total stage steps)
%
% Outputs: t = tspan
%          y = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%              y(t*) is a column vector of length m.
%
% Note: to run in fixed-step mode, just call the solver using hmin=hmax as
% the desired time step size.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

% extract ERK method information from B
[Brows, Bcols] = size(B);
s = Bcols - 1;
c = B(1:s,1);
b = (B(s+1,2:s+1))';
A = B(1:s,2:s+1);
if (Brows > Bcols)
   embedded = 1;
   b2 = (B(s+2,2:s+1))';
else
   embedded = 0;
end

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% initialize diagnostics
h_a = 0;
h_s = 0;

% set the solver parameters
dt_safety  = 0.1;
dt_reduce  = 0.1;
dt_stable  = 0.5;


% store temporary states
t = tvals(1);
Ynew = Y0;

% create Fdata structure (since used in nonlinear solver, set IRK values)
Fdata.fname = fcn;
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
      
      % bound internal time step using inputs
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time
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
	 
	 % construct stage solution
	 %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
	 z(:,stage) = Y0;
	 for j=1:stage-1
	    z(:,stage) = z(:,stage) + h*A(stage,j)*feval(fcn,t+h*c(j),z(:,j));
	 end
	 nsteps = nsteps + 1;
	 
      end

      % use stage solutions to determine time-evolved answer, update old
      Fdata.t = t;
      [Ynew,Y2] = Y_ERK(z,Fdata);
      Y0 = Ynew;
      
      % update time, work counter
      t = t + h;
      
      % for embedded methods, estimate error and update time step,
      % assuming that method local truncation error order equals number of stages
      if (embedded) 
	 h = h_estimate(Ynew, Y2, h, rtol, atol, s);
      else
	 h = hmin;
      end
      
      % limit time step by stability condition
      hstab = dt_stable * feval(EStabFn, t, Ynew);
      if (h < hstab)
	 h_a = h_a + s;
      else
	 h_s = h_s + s;
	 %fprintf('   h limited for explicit stability: h_i = %g, h_e = %g\n',h,hstab);
      end
      h = min([h, hstab]);
	 
   end

   % store updated solution
   Y(:,tstep) = Ynew;
   
end


% output diagnostics
%fprintf('solve_ERK: took %i accuracy and %i stability steps\n',h_a,h_s);


% end function
