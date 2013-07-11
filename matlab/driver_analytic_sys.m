% driver for stiff ODE system with analytical solution:
%    u' = A*u
% where u(0) = [1,1,1]', A = V*D*Vinv, 
%      V = [1 -1 1; -1 2 1; 0 -1 2];
%      Vinv = 0.25*[5 1 -3; 2 2 -2; 1 1 1];
%      D = [-0.5 0 0; 0 -0.1 0; 0 0 lam];
% where lam is a large negative number. The analytical solution to
% this problem is 
%   Y(t) = V*exp(D*t)*Vinv*Y0
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

clear

% set problem parameters
lam = -100;

% set the total integration time
Tf = 0.05;

% set desired output times
tout = linspace(0,Tf,11);

% set the time step size bounds, tolerance
hmin = 1e-12;
hmax = 1.0;
rtol = 1e-6;
atol = 1e-10*ones(1,1);
hmethod = 3;

% get the DIRK Butcher tables
mname = 'SDIRK-5-4';
B = butcher(mname);


% set the initial conditions
Y0 = [1,1,1]';

% store problem parameters
global Pdata;
Pdata.lam = lam;


fprintf('\nAnalytical test (rtol = %g, atol = %g)\n',rtol,atol(1))

% get "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s('f_analytic_sys', tout, Y0, opts);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK('f_analytic_sys', 'J_analytic_sys', tout, Y0, B, rtol, ...
    atol, hmin, hmax, hmethod);

% compute error
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sum(sum((Y'-Ytrue).^2));
err_rms = sqrt(err_rms/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);



% end of script
