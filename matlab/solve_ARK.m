function [tvals,Y,nsteps] = solve_ARK(fcnI, fcnE, Jfcn, EStabFn, tvals, Y0, Bi, Be, tol, hmin, hmax)
% usage: [tvals,Y,nsteps] = solve_ARK(fcnI, fcnE, Jfcn, EStabFn, tvals, Y0, Bi, Be, tol, hmin, hmax)
%
% Adaptive time step additive RK solver for the vector-valued ODE problem
%     y' = F(t,Y), t in tspan,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:  fcnI = function name for implicit portion of ODE right-hand side
%          fcnE = function name for explicit portion of ODE right-hand side,
%                 where   F(t,Y) = fcnI(t,Y) + fcnE(t,Y)
%          Jfcn = function name for Jacobian of fcnI, J(t,Y)
%          EStabFn = function name for stability constraint on fcnE
%          tvals = [t0, t1, t2, ..., tN]
%          Y0 = initial values
%          Bi = Butcher matrix for SDIRK coefficients, of the form
%                Bi = [c A;
%                      0 b;
%                      0 b2 ]
%          Be = Butcher matrix for SDIRK coefficients, of the form
%                Be = [c A;
%                      0 b;
%                      0 b2 ]
%              The b2 row is optional in both Bi and Be, but if not provided
%              the method will default to taking fixed step sizes, at
%              min(hmin, hstab), where hstab is the result of EStabFn(t,Y).
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

% extract DIRK method information from Bi
[Brows, Bcols] = size(Bi);
si = Bcols - 1;
ci = Bi(1:si,1);
bi = (Bi(si+1,2:si+1))';
Ai = Bi(1:si,2:si+1);
if (Brows > Bcols) 
   embedded = 1;
   bi2 = (Bi(si+2,2:si+1))';
else
   embedded = 0;
end

% extract ERK method information from Be
[Brows, Bcols2] = size(Be);
se = Bcols2 - 1;
ce = Be(1:se,1);
be = (Be(se+1,2:se+1))';
Ae = Be(1:se,2:se+1);
if ((Brows > Bcols2) && (embedded == 1))
   embedded = 1;
   be2 = (Be(si+2,2:si+1))';
else
   embedded = 0;
end

% ensure that DIRK and ERK method match at internal stage times
if (si == se-1)
   % if DIRK method not written with initial explicit stage, pad DIRK with a
   % zero first stage to match ERK
   Ai = [0, 0*bi; 0*ci; Ai];
   ci = [0; ci];
   bi = [0, bi];
   si = si + 1;
   if (embedded)
      bi2 = [0, bi2];
      Bi = [ci, Ai; 0, bi; 0, bi2];
   else
      Bi = [ci, Ai; 0, bi];
   end
end

% ensure that ci now equals ce, otherwise return with error
if (si ~= se)
   error('incompatible SDIRK/EKR pair for ARK method, internal stage mismatch!')
end
if (norm(ci-ce) > sqrt(eps))
   error('incompatible SDIRK/EKR pair for ARK method, internal stage mismatch!')
end


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
dt_stable  = 0.5;


% store temporary states
t = tvals(1);
Ynew = Y0;

% create Fdata structure (since used in nonlinear solver, set IRK values)
Fdata.fname = fcnI;
Fdata.Jname = Jfcn;
Fdata.B = Bi;
Fdata.s = si;

% add in explicit method data for eventual solution
Fdata.fnameE = fcnE;
Fdata.Be = Be;

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
      z = zeros(m,si);

      % reset stage success flag
      st_success = 0;
      
      % loop over stages
      for stage=1:si
	 
	 % set stage initial guess as previous stage solution
	 Yguess = Ynew;
      
	 % set stage number into Fdata structure
	 Fdata.stage = stage;
	 
	 % construct RHS comprised of old time data
	 %    zi = y_n + h*sum_{j=1}^i (Ai(i,j)*fi(zj)) 
	 %             + h*sum_{j=1}^{i-1} (Ae(i,j)*fe(zj))
	 %    zi = y_n + h*Ai(i,i)*fi(zi)
	 %             + h*sum_{j=1}^{i-1} (Ai(i,j)*fi(zj) + Ae(i,j)*fe(zj))
	 %    zi - h*Ai(i,i)*fi(zi) = 
	 %             y_n + h*sum_{j=1}^{i-1} (Ai(i,j)*fi(zj) + Ae(i,j)*fe(zj))
	 % =>
	 %    rhs = y_n + h*sum_{j=1}^{i-1} (Ai(i,j)*fi(zj) + Ae(i,j)*fe(zj))
	 Fdata.rhs = Y0;
	 for j=1:stage-1
	    Fdata.rhs = Fdata.rhs + h*Ai(stage,j)*feval(fcnI,t+h*ci(j),z(:,j)) ...
		                  + h*Ae(stage,j)*feval(fcnE,t+h*ce(j),z(:,j));
	 end
	 
	 % call newton solver to compute new stage solution
	 Fdata.t = t;
	 [Ynew,ierr] = newton_damped('F_DIRK', 'A_DIRK', Yguess, Fdata, ...
	     newt_tol, newt_maxit, newt_alpha);
	 nsteps = nsteps + 1;
	 
	 % check newton error flag, if failure break out of stage loop
	 if (ierr ~= 0) 
	    fprintf('solve_ARK warning: stage failure, cutting timestep\n');
	    st_success = 1;
	    break;
	 end
	 
	 % update stored solution with new value
	 z(:,stage) = Ynew;
	 
      end

      % check whether all stages solved successfully
      if (st_success == 0) 
      
	 % use stage solutions to determine time-evolved answer, update old
	 [Ynew,Y2] = Y_ARK(z,Fdata);
	 Y0 = Ynew;

	 % update time, work counter
	 t = t + h;
   
	 % for embedded methods, estimate error and update time step,
	 % assuming that method order equals number of stages-1 (this should
	 % be an input argument) 
	 if (embedded) 
	    err = norm(Y2-Ynew,inf) + sqrt(eps)*tol;
	    h = max([h * (dt_safety*tol/err)^(1/(si-1)), hmin]);
	 else
	    h = hmin;
	 end
	 
	 % limit time step by explicit stability condition
	 hstab = dt_stable * feval(EStabFn, t, Ynew);
	 h = min([h, hstab]);
	 
      % if any stages failed, reduce step size and retry
      else
	 
	 % reset solution guess, update work counter, reduce time step
	 Ynew = Y0;
	 h = h * dt_reduce;
      
      end
      
   end

   % store updated solution
   Y(:,tstep) = Ynew;
   
end


% end function
