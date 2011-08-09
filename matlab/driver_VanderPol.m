% driver for stiff ODE version of Van der Pol test problem:
%    u' = v
%    v' = (v - v*u^2 - u)/ep
% where u(0) = 2,  v(0) = -0.6666654321121172, and ep = 1e-5.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear

% set problem parameters
ep = 1e-4;

% set the total integration time
Tf = 1.5;

% set desired output times
tout = [0,Tf];

% set the time step size bounds, tolerance
%hmin = 10^(-4.5);
hmin = 1e-6;
hmax = 0.1;
tol = 1e-6;

% get the DIRK Butcher table
% mname = 'ARK3(2)4L[2]SA-ESDIRK';
mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);

% set the initial conditions
u0 = -2;
v0 = -2.3553013976081;
Y0 = [u0; v0];

% store problem parameters
global Pdata;
Pdata.ep = ep;

%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning van der Pol test with integrator: %s  (tol = %g)\n',mname,tol)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK2('f_VanderPol', 'J_VanderPol', tout, Y0, B, tol, hmin, hmax);

% get "true" solution
fprintf('\nComputing "true" solution with ode15s\n')
opts = odeset('RelTol',1e-10, 'AbsTol',1e-14*ones(size(Y0)),...
              'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_VanderPol', tout, Y0, opts);

% compute error
err_max = 0;
err_rms = 0;
for j=1:length(Y0)
   diff = (Y(j,end) - Ytrue(end,j))/Ytrue(end,j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/length(Y0));
fprintf('\nAccuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);


% end of script
