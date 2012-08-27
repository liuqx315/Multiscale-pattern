% driver for stiff ODE version of analytical test problem:
%    u' = lam*u + 1/(1+t^2) - lam*atan(t)
% where u(0) = 0, lam = -1/ep, and ep = 1e-2.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% May 2012
% All Rights Reserved

clear

% set problem parameters
ep = 1e-2;

% set the total integration time
Tf = 10.0;

% set desired output times
tout = linspace(0,Tf,11);


% set the time step size bounds, tolerance
hmin = 1e-12;
hmax = 1.0;
rtol = 1e-6;
atol = 1e-10*ones(1,1);
hmethod = 1;

% get the DIRK Butcher tables
mname = 'Billington-SDIRK';
B = butcher(mname);


% set the initial conditions
Y0 = 0;

% store problem parameters
global Pdata;
Pdata.ep = ep;


fprintf('\nAnalytical test (rtol = %g, atol = %g)\n',rtol,atol(1))

% get true solution
Ytrue = atan(tout);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK('f_analytic', 'J_analytic', tout, Y0, B, rtol, ...
    atol, hmin, hmax, hmethod);

% compute error
err_max = max(max(abs(Y-Ytrue)));
err_rms = sum(sum((Y-Ytrue).^2));
err_rms = sqrt(err_rms/(numel(Y)-1));  % subtract 1 since IC is exact
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);



% end of script
