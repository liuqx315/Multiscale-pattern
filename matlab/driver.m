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

% set the time step size
h = 0.1;

% $$$ % set the IRK table of coefficients (3 stage Gauss method of order 6)
% $$$ mname = 'Gauss 3-stage';
% $$$ c = [0.5-sqrt(15)/10; 0.5; 0.5+sqrt(15)/10];
% $$$ b = [5; 8; 5]./18;
% $$$ A = [5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30; ...
% $$$      5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24; ...
% $$$      5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];

% set the IRK table of coefficients (5 stage Radau IIA method of order 9)
mname = 'RadauIIA 5-stage';
c = [0.05710419611451768; 0.2768430136381238; 0.5835904323689168; ...
     0.8602401356562195; 1];
b = [0.1437135607912259; 0.2813560151494621; 0.3118265229757413; ...
     0.2231039010835707; 0.04]; 
A = [0.07299886431790337, -0.02673533110794565, 0.01867692976398445, ...
       -0.01287910609330652, 0.005042839233882052; 
     0.1537752314791824, 0.1462148678474935, -0.03644456890512816, ...
0.02123306311930480, -0.007935579902728813;
     0.1400630456848099, 0.2989671294912833, 0.1675850701352492, ...
       -0.03396910168661794, 0.01094428874419233; 
     0.1448943081095342, 0.2765000687601608, 0.3257979229104191, ...
0.1287567532549115, -0.01570891737880607; 
     0.1437135607912259, 0.2813560151494621, 0.3118265229757413, ...
        0.2231039010835707, 0.04];

% store the Butcher table in a single matrix
B = [c, A; 0, b'];


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
[t,Y,ns] = solve_IRK('fbruss1', 'Jbruss1', tout, Y0, B, h);
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
[t,Y,ns] = solve_IRK('fbruss2', 'Jbruss2', tout, Y0, B, h);
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
[t,Y,ns] = solve_IRK('fbruss3', 'Jbruss3', tout, Y0, B, h);
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


%%%%%% Second set of tests -- speed %%%%%%%

%%%%%%%%%%%%%%%%
% Test 1: initially very fast, then slower evolution
% a = 1.2; 
% b = 2.5; 
% ep = 1e-5;
u0 = 3.9;
v0 = 1.1;
w0 = 2.8;
Y0 = [u0; v0; w0];

% integrate using student's solver for various Tf
ts = cputime;
tout = [0,1];
[t,Y,ns] = solve_IRK('fbruss1', 'Jbruss1', tout, Y0, B, h);
tout = [0,10];
[t,Y,ns2] = solve_IRK('fbruss1', 'Jbruss1', tout, Y0, B, h);
tout = [0,50];
[t,Y,ns3] = solve_IRK('fbruss1', 'Jbruss1', tout, Y0, B, h);
tf = cputime;
fprintf('Speed Test 1 time = %.5e,  work = %i\n',tf-ts,ns+ns2+ns3)
fprintf('\n');


%%%%%%%%%%%%%%%%
% Test 2: initially slow then fast evolution
% a = 1; 
% b = 3.5; 
% ep = 5e-6;
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];

% integrate using student's solver for various Tf
ts = cputime;
tout = [0,1];
[t,Y,ns] = solve_IRK('fbruss2', 'Jbruss2', tout, Y0, B, h);
tout = [0,10];
[t,Y,ns2] = solve_IRK('fbruss2', 'Jbruss2', tout, Y0, B, h);
tout = [0,50];
[t,Y,ns3] = solve_IRK('fbruss2', 'Jbruss2', tout, Y0, B, h);
tf = cputime;
fprintf('Speed Test 2 time = %.5e,  work = %i\n',tf-ts,ns+ns2+ns3)
fprintf('\n');


%%%%%%%%%%%%%%%%
% Test 3: initially very rapid, then equilibrium or steady growth
% a = 0.5; 
% b = 3; 
% ep = 5e-4;
u0 = 3;
v0 = 3;
w0 = 3.5;
Y0 = [u0; v0; w0];

% integrate using student's solver for various Tf
ts = cputime;
tout = [0,1];
[t,Y,ns] = solve_IRK('fbruss3', 'Jbruss3', tout, Y0, B, h);
tout = [0,10];
[t,Y,ns2] = solve_IRK('fbruss3', 'Jbruss3', tout, Y0, B, h);
tout = [0,50];
[t,Y,ns3] = solve_IRK('fbruss3', 'Jbruss3', tout, Y0, B, h);
tf = cputime;
fprintf('Speed Test 3 time = %.5e,  work = %i\n',tf-ts,ns+ns2+ns3)
fprintf('\n');


% end of script
