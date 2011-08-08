% driver for 1D Cahn-Hilliard problem,
%   u_t = -\partial_xx (c^2 \partial_xx u - u(u^2 - 1)), on [0,1]
%      u_x = 0     at x=0,x=1
%      u_xxx = 0   at x=0,x=1
%      u(x,0) = u_0 + r,  r is a uniform random number in [-0.05,0.05].
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear
%maxNumCompThreads(1);

% set the problem parameters
c = 0.00078;

% set the total integration time, mesh size
%Tf = 0.1;
Tf = 0.001;
N = 2049;

% set desired output times
tout = linspace(0,Tf,50);

% set the time step size bounds, tolerance
hmin = 1e-12;
hmax = 1e-5;
tol  = 1e-6;

% get the DIRK Butcher table
mname = 'ARK3(2)4L[2]SA-ESDIRK';
%mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);

% store the Butcher table in a single matrix
B = [c, A; 0, b'; 0, b2'];

% initial conditions
%      u(x,0) = u_0 + r,  r is a uniform random number in [-0.05,0.05].
Y0 = 0.5 + 0.01*(rand(N,1)-0.5);

% boundary conditions
%      u_x = 0     at x=0,x=1
%      u_xxx = 0   at x=0,x=1
Y0(1)   = 39/17*Y0(3)   - 28/17*Y0(4)   + 6/17*Y0(5);
Y0(2)   = 67/34*Y0(3)   - 21/17*Y0(4)   + 9/34*Y0(5);
Y0(N-1) = 67/34*Y0(N-2) - 21/17*Y0(N-3) + 9/34*Y0(N-4);
Y0(N)   = 39/17*Y0(N-2) - 28/17*Y0(N-3) + 6/17*Y0(N-4);

% store the problem parameters for function evaluations
Pdata.c = c;
Pdata.dx = 1/(N-1);
Pdata.n = N;

printf('\nRunning 1D Cahn-Hilliard test with integrator: ')
disp(mname)
disp('   ')

% integrate using built-in method and tight tolerances for "true" solution
Ytrue = Y0;

% integrate using solver, and compare with "true" solution
err_max = 0;
err_rms = 0;
[t,Y,ns] = solve_DIRK2('fCH_1D', 'JCH_1D', tout, Y0, B, tol, hmin, hmax);
for j=1:2*N
   diff = (Y(j,end) - Ytrue(j))/Ytrue(j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/(2*N));
fprintf('Accuracy Test Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);



% end of script
