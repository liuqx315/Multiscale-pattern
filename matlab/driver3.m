% driver for racing part of projective integration project
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2011
% All Rights Reserved

clear
%maxNumCompThreads(1);

% set the total integration time
Tf = 10;

% set desired output times
tout = [0,Tf];

% set the time step size bounds, tolerance
hmin = 1e-6;
hmax = 1.0;
tol = 1e-6;

% $$$ % set the DIRK table of coefficients -- ARK3(2)4L[2]SA-ESDIRK
% $$$ mname = 'ARK3(2)4L[2]SA-ESDIRK';
% $$$ c = [0; 1767732205903/2027836641118; 3/5; 1];
% $$$ b = [1471266399579/7840856788654; -4482444167858/7529755066697; 
% $$$      11266239266428/11593286722821; 1767732205903/4055673282236];
% $$$ b2 = [2756255671327/12835298489170; -10771552573575/22201958757719; 
% $$$       9247589265047/10645013368117; 2193209047091/5459859503100];
% $$$ gamma = 1767732205903/4055673282236;
% $$$ A = [0, 0, 0, 0;
% $$$      1767732205903/4055673282236, gamma, 0, 0;
% $$$      2746238789719/10658868560708, -640167445237/6845629431997, gamma, 0;
% $$$      1471266399579/7840856788654, -4482444167858/7529755066697, 11266239266428/11593286722821, gamma];

% set the DIRK table of coefficients -- ARK4(3)6L[2]SA-ESDIRK
mname = 'ARK4(3)6L[2]SA-ESDIRK';
c = [0; 1/2; 83/250; 31/50; 17/20; 1];
b = [82889/524892; 0; 15625/83664; 69875/102672; -2260/8211; 1/4];
b2 = [4586570599/29645900160; 0; 178811875/945068544; 814220225/1159782912; ...
   -3700637/11593932; 61727/225920];
gamma = 1/4;
A = [0, 0, 0, 0, 0, 0;
     1/4, gamma, 0, 0, 0, 0;
     8611/62500, -1743/31250, gamma, 0, 0, 0;
     5012029/34652500, -654441/2922500, 174375/388108, gamma, 0, 0;
     15267082809/155376265600, -71443401/120774400, 730878875/902184768, 2285395/8070912, gamma, 0;
     82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4];

% store the Butcher table in a single matrix
B = [c, A; 0, b'; 0, b2'];


printf('\nRunning tests with integrator: ')
disp(mname)
disp('   ')

%%%%%% First set of tests -- accuracy %%%%%%%

%%%%%%%%%%%%%%%%
% Test 1: initially very fast, then slower evolution
% a = 1.2; 
% b = 2.5; 
% ep = 1e-5;
u0 = 3.9;
v0 = 1.1;
w0 = 2.8;
Y0 = [u0; v0; w0];

% integrate using student's solver, and compare with "true" solution
Ytrue = load('Ytrue_test1');
err_max = 0;
err_rms = 0;
[t,Y,ns] = solve_DIRK2('fbruss1', 'Jbruss1', tout, Y0, B, tol, hmin, hmax);
for j=1:3
   diff = (Y(j,end) - Ytrue(j))/Ytrue(j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/3);
fprintf('Accuracy Test 1 Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);


%%%%%%%%%%%%%%%%
% Test 2: initially slow then fast evolution
% a = 1; 
% b = 3.5; 
% ep = 5e-6;
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];

% integrate using student's solver, and compare with "true" solution
Ytrue = load('Ytrue_test2');
err_max = 0;
err_rms = 0;
[t,Y,ns] = solve_DIRK2('fbruss2', 'Jbruss2', tout, Y0, B, tol, hmin, hmax);
for j=1:3
   diff = (Y(j,end) - Ytrue(j))/Ytrue(j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/3);
fprintf('Accuracy Test 2 Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);


%%%%%%%%%%%%%%%%
% Test 3: initially very rapid, then equilibrium or steady growth
% a = 0.5; 
% b = 3; 
% ep = 5e-4;
u0 = 3;
v0 = 3;
w0 = 3.5;
Y0 = [u0; v0; w0];

% integrate using student's solver, and compare with "true" solution
Ytrue = load('Ytrue_test3');
err_max = 0;
err_rms = 0;
[t,Y,ns] = solve_DIRK2('fbruss3', 'Jbruss3', tout, Y0, B, tol, hmin, hmax);
for j=1:3
   diff = (Y(j,end) - Ytrue(j))/Ytrue(j);
   err_max = max([err_max, abs(diff)]);
   err_rms = err_rms + diff^2;
end
err_rms = sqrt(err_rms/3);
fprintf('Accuracy Test 3 Results:\n')
fprintf('   maxerr = %.5e,   rmserr = %.5e\n',err_max,err_rms);
fprintf('   work = %i\n',ns);
fprintf('\n');



% end of script
