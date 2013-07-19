function butcher_test_method(B,mname)
% function to check accuracy of a given RK Butcher table
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% May 2012
% All Rights Reserved

% set the total integration time
Tf = 10;

% set desired output times
tout = [0,Tf];

% set the time step sizes
hvals = logspace(-1,-3,10);

% initial guess
u0 = 1;
v0 = 0;
w0 = 0.9;
Y0 = [u0; v0; w0];

% get "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',1e-14,'InitialStep',1e-12, 'MaxStep',1e-3);
[t,Ytrue] = ode45('f_test', tout, Y0, opts);

% determine method type (ERK,DIRK,IRK)
[Brows,Bcols] = size(B);
dirk = 0;
irk = 0;
s = Bcols-1;
A = B(1:s,2:Bcols);
for i = 1:s
   if (abs(A(i,i)) > 1e-10)
      dirk = 1;
   end
   if (i < s)
      if (sum(abs(A(i,i+1:s))) > 1e-10)
         irk = 1;
      end
   end
end

% initialize result array
m = zeros(length(hvals),2);

% iterate over h values
for i = 1:length(hvals)
   h = hvals(i);
   if (irk == 1) 
      [t,Y,ns] = solve_IRK('f_test','J_test',tout,Y0,B,1e2,1e2,h,h,1);
   elseif (dirk == 1) 
      [t,Y,ns] = solve_DIRK('f_test','J_test',tout,Y0,B,1e2,1e2,h,h,1);
   else
      [t,Y,ns] = solve_ERK('f_test','EStab_test',tout,Y0,B,1e2,1e2,h,h,1);
   end
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(i,1:2) = [ norm(diff,inf) , norm(diff) ];
end

% output results for method
fprintf('\n  Results for %s:\n',mname);
fprintf('      h           errI        err2\n');
fprintf('  --------------------------------------\n')
for ih = 1:length(hvals)
   errI(ih) = m(ih,1);
   err2(ih) = m(ih,2);
   fprintf('  %10.2e  %10.2e  %10.2e\n',hvals(ih),errI(ih),err2(ih))
end
fprintf('  --------------------------------------\n')
p1 = polyfit(log(hvals(2:end-2)),log(errI(2:end-2)),1);
p2 = polyfit(log(hvals(2:end-2)),log(err2(2:end-2)),1);
p = (p1(1)+p2(1))/2;
fprintf('  Measured order of accuracy = %g\n',p);



%--------------------
% If an embedding is available, test that as well
if (Brows == Bcols)
   return
end

% swap b arrays
btmp = B(Brows,:);
B(Brows,:) = B(Brows-1,:);
B(Brows-1,:) = btmp;

% perform test again
m = zeros(length(hvals),2);
for i = 1:length(hvals)
   h = hvals(i);
   if (irk == 1) 
      [t,Y,ns] = solve_IRK('f_test','J_test',tout,Y0,B,1e-2,1e-2,h,h,1);
   elseif (dirk == 1) 
      [t,Y,ns] = solve_DIRK('f_test','J_test',tout,Y0,B,1e-2,1e-2,h,h,1);
   else
      [t,Y,ns] = solve_ERK('f_test','EStab_test',tout,Y0,B,1e-2,1e-2,h,h,1);
   end
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(i,1:2) = [ norm(diff,inf) , norm(diff) ];
end
fprintf('\n  Results for %s embedding:\n',mname);
fprintf('      h           errI        err2\n');
fprintf('  --------------------------------------\n')
for ih = 1:length(hvals)
   errI(ih) = m(ih,1);
   err2(ih) = m(ih,2);
   fprintf('  %10.2e  %10.2e  %10.2e\n',hvals(ih),errI(ih),err2(ih))
end
fprintf('  --------------------------------------\n')
p1 = polyfit(log(hvals(2:end-2)),log(errI(2:end-2)),1);
p2 = polyfit(log(hvals(2:end-2)),log(err2(2:end-2)),1);
p = (p1(1)+p2(1))/2;
fprintf('  Measured order of accuracy = %g\n',p);



% end of function
