% driver for 2D radiative ionization PDE test
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear
%maxNumCompThreads(1);

% set the total integration time, mesh size
Tf = 80;
N = 1.0/0.002;

% set the brusselator parameters
alpha = 0.6; 
beta = 2; 

% set desired output times
tout = linspace(0,Tf,50);

% set the time step size bounds, tolerance
hmin = 1e-6;
hmax = 1.0;
tol = 1e-6;

% get the DIRK Butcher table
% mname = 'ARK3(2)4L[2]SA-ESDIRK';
mname = 'ARK4(3)6L[2]SA-ESDIRK';
B = butcher(mname);


printf('\nRunning 1D brusselator tests with integrator: ')
disp(mname)
disp('   ')

% initial conditions
xspan = linspace(0,1,N)';
u0 = alpha + 0.1*sin(pi*xspan);
v0 = beta/alpha + 0.1*sin(pi*xspan);
Y0 = [u0; v0];

% set problem parameter data structure
Pdata.N = N;

% integrate using built-in method and tight tolerances for "true" solution


% integrate using solver, and compare with "true" solution
err_max = 0;
err_rms = 0;
[t,Y,ns] = solve_DIRK2('fbruss_1D', 'Jbruss_1D', Pdata, tout, Y0, B, tol, hmin, hmax);
%for j=1:2*N
%   diff = (Y(j,end) - Ytrue(j))/Ytrue(j);
%   err_max = max([err_max, abs(diff)]);
%   err_rms = err_rms + diff^2;
%end
%err_rms = sqrt(err_rms/(2*N));
fprintf('Accuracy Test Results:\n')
%fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);



% end of script
