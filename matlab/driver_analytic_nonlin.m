% driver for stiff ODE version of analytical test problem:
%    y' = (t+1)*exp(-y)
% where y(0) = 0.  This has analytical solution y(t) = log(0.5*((t+1)^2+1))
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

clear

% set the total integration time
Tf = 10;

% set desired output times
tout = linspace(0,Tf,11);

% set the time step size bounds, tolerance
%hmin = 1e-12;
%hmax = 1.0;
%rtol = 1e-6;
%atol = 1e-10*ones(1,1);
hmin = 1.0;
hmax = 1.0;
rtol = 1.0;
atol = 1.0*ones(1,1);
hmethod = 3;

% get the DIRK Butcher tables
mname = 'SDIRK-5-4';
B = butcher(mname);

% set the initial conditions
Y0 = 0;


fprintf('\nNonlinear analytical test (rtol = %g, atol = %g)\n',rtol,atol(1))

% get true solution
Ytrue = log(((tout+1).^2+1)/2);


%%%%%%%%%%%%%%%%%%%%
fprintf('\nRunning with SDIRK integrator: %s\n',mname)

% integrate using adaptive solver
[t,Y,ns] = solve_DIRK('f_analytic_nonlin', 'J_analytic_nonlin', ...
    tout, Y0, B, rtol, atol, hmin, hmax, hmethod);

% compute error
err_max = max(max(abs(Y-Ytrue)));
err_rms = sum(sum((Y-Ytrue).^2));
err_rms = sqrt(err_rms/numel(Y));
fprintf('\nResults:\n')
fprintf('        t           y         error\n');
fprintf('   ------------------------------------\n');
for i=1:length(Y)
   fprintf('  %10.6f  %10.6f  %12.5e\n', tout(i), Y(i), Ytrue(i)-Y(i));
end
fprintf('   ------------------------------------\n');

fprintf('\nOverall Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);



% end of script
