% driver to check accuracy of various RK methods in butcher.m
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2011
% All Rights Reserved

clear
!\rm butcher_tests.txt
diary butcher_tests.txt

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
fprintf('Calling ode45 for reference solution\n');
opts = odeset('RelTol',1e-12, 'AbsTol',1e-14,'InitialStep',1e-12, 'MaxStep',1e-3);
[t,Ytrue] = ode45('f_test', tout, Y0, opts);
figure(1)
plot(t,Ytrue(:,1),t,Ytrue(:,2),t,Ytrue(:,3))
xlabel('time')
legend('u','v','w')
title('Reference Solution')

% iterate over h values
fprintf('\nRunning tests with h = ');
for i = 1:length(hvals)
   
   % set h value, refresh method integer
   h = hvals(i);
   imethod = 0;
   
   fprintf(' %g, ',h);

   %%%%%%%%%%%%%%%%%%% ERK Methods %%%%%%%%%%%%%%%%%%%

   %    ARK3(2)4L[2]SA-ERK
   imethod = imethod + 1;
   mname = 'ARK3(2)4L[2]SA-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  ARK3(2)E' , norm(diff,inf) , norm(diff) };

   %    ARK4(3)6L[2]SA-ERK
   imethod = imethod + 1;
   mname = 'ARK4(3)6L[2]SA-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  ARK4(3)E' , norm(diff,inf) , norm(diff) };
   
   %    ARK5(4)8L[2]SA-ERK
   imethod = imethod + 1;
   mname = 'ARK5(4)8L[2]SA-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  ARK5(4)E' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(2,3,3)-ERK
   imethod = imethod + 1;
   mname = 'Ascher(2,3,3)-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(2,3,3)E' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(2,3,2)-ERK
   imethod = imethod + 1;
   mname = 'Ascher(2,3,2)-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(2,3,2)E' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(2,2,2)-ERK
   imethod = imethod + 1;
   mname = 'Ascher(2,2,2)-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(2,2,2)E' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(3,4,3)-ERK
   imethod = imethod + 1;
   mname = 'Ascher(3,4,3)-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(3,4,3)E' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(4,4,3)-ERK
   imethod = imethod + 1;
   mname = 'Ascher(4,4,3)-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(4,4,3)E' , norm(diff,inf) , norm(diff) };

   %    Cooper4-ERK
   imethod = imethod + 1;
   mname = 'Cooper4-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  Cooper4E' , norm(diff,inf) , norm(diff) };
   
   %    Cooper6-ERK
   imethod = imethod + 1;
   mname = 'Cooper6-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  Cooper6E' , norm(diff,inf) , norm(diff) };

   %    Heun-Euler-ERK
   imethod = imethod + 1;
   mname = 'Heun-Euler-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { 'Heun-Euler' , norm(diff,inf) , norm(diff) };

   %    Bogacki-Shampine-ERK
   imethod = imethod + 1;
   mname = 'Bogacki-Shampine-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  BogShamp' , norm(diff,inf) , norm(diff) };

   %    Fehlberg-ERK
   imethod = imethod + 1;
   mname = 'Fehlberg-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  Fehlberg' , norm(diff,inf) , norm(diff) };

   %    Cash-Karp-ERK
   imethod = imethod + 1;
   mname = 'Cash-Karp-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  CashKarp' , norm(diff,inf) , norm(diff) };

   %    Dormand-Prince-ERK
   imethod = imethod + 1;
   mname = 'Dormand-Prince-ERK';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { 'DormPrince' , norm(diff,inf) , norm(diff) };

   %    ERK-1-1
   imethod = imethod + 1;
   mname = 'ERK-1-1';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '   ERK-1-1' , norm(diff,inf) , norm(diff) };

   %    ERK-2-2
   imethod = imethod + 1;
   mname = 'ERK-2-2';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '   ERK-2-2' , norm(diff,inf) , norm(diff) };

   %    ERK-3-3
   imethod = imethod + 1;
   mname = 'ERK-3-3';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '   ERK-3-3' , norm(diff,inf) , norm(diff) };

   %    ERK-4-4
   imethod = imethod + 1;
   mname = 'ERK-4-4';
   B = butcher(mname);
   [t,Y,ns] = solve_ERK('f_test', 'EStab_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '   ERK-4-4' , norm(diff,inf) , norm(diff) };

   
   %%%%%%%%%%%%%%%%%%% DIRK Methods %%%%%%%%%%%%%%%%%%%

   
   %    ARK3(2)4L[2]SA-ESDIRK
   imethod = imethod + 1;
   mname = 'ARK3(2)4L[2]SA-ESDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  ARK3(2)I' , norm(diff,inf) , norm(diff) };

   %    ARK4(3)6L[2]SA-ESDIRK
   imethod = imethod + 1;
   mname = 'ARK4(3)6L[2]SA-ESDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  ARK4(3)I' , norm(diff,inf) , norm(diff) };
   
   %    ARK5(4)8L[2]SA-ESDIRK
   imethod = imethod + 1;
   mname = 'ARK5(4)8L[2]SA-ESDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  ARK5(4)I' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(2,3,3)-SDIRK
   imethod = imethod + 1;
   mname = 'Ascher(2,3,3)-SDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(2,3,3)I' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(3,4,3)-SDIRK
   imethod = imethod + 1;
   mname = 'Ascher(3,4,3)-SDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(3,4,3)I' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(2,3,2)-SDIRK
   imethod = imethod + 1;
   mname = 'Ascher(2,3,2)-SDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(2,3,2)I' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(2,2,2)-SDIRK
   imethod = imethod + 1;
   mname = 'Ascher(2,2,2)-SDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(2,2,2)I' , norm(diff,inf) , norm(diff) };
   
   %    Ascher(4,4,3)-SDIRK
   imethod = imethod + 1;
   mname = 'Ascher(4,4,3)-SDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' A(4,4,3)I' , norm(diff,inf) , norm(diff) };
   
   %    Cooper4-ESDIRK
   imethod = imethod + 1;
   mname = 'Cooper4-ESDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  Cooper4I' , norm(diff,inf) , norm(diff) };
   
   %    Cooper6-ESDIRK
   imethod = imethod + 1;
   mname = 'Cooper6-ESDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  Cooper6I' , norm(diff,inf) , norm(diff) };
   
   %    TRBDF2-ESDIRK
   imethod = imethod + 1;
   mname = 'TRBDF2-ESDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '    TRBDF2' , norm(diff,inf) , norm(diff) };
   
   %    TRX2-ESDIRK
   imethod = imethod + 1;
   mname = 'TRX2-ESDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '      TRX2' , norm(diff,inf) , norm(diff) };
   
   %    Billington-SDIRK
   imethod = imethod + 1;
   mname = 'Billington-SDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { 'Billington' , norm(diff,inf) , norm(diff) };

   %    Cash(5,2,4)-SDIRK
   imethod = imethod + 1;
   mname = 'Cash(5,2,4)-SDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  C(5,2,4)' , norm(diff,inf) , norm(diff) };
   
   %    Cash(5,3,4)-SDIRK
   imethod = imethod + 1;
   mname = 'Cash(5,3,4)-SDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  C(5,3,4)' , norm(diff,inf) , norm(diff) };
   
   %    Ismail(7,4,5)-ESDIRK
   imethod = imethod + 1;
   mname = 'Ismail(7,4,5)-ESDIRK';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  I(7,4,5)' , norm(diff,inf) , norm(diff) };
   
   %    SDIRK-2-2
   imethod = imethod + 1;
   mname = 'SDIRK-2-2';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' SDIRK 2-2' , norm(diff,inf) , norm(diff) };
   
   %    SDIRK-4-5
   imethod = imethod + 1;
   mname = 'SDIRK-4-5';
   B = butcher(mname);
   [t,Y,ns] = solve_DIRK('f_test', 'J_test', tout, Y0, B, 1e-2, 1e-2, h, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' SDIRK 4-5' , norm(diff,inf) , norm(diff) };

   %%%%%%%%%%%%%%%%%%% IRK Methods %%%%%%%%%%%%%%%%%%%

   %    IRK-1-1
   imethod = imethod + 1;
   mname = 'IRK-1-1';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '   IRK 1-1' , norm(diff,inf) , norm(diff) };

   %    Crank-Nicolson-2-2-IRK
   imethod = imethod + 1;
   mname = 'Crank-Nicolson-2-2-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '    CN 2-2' , norm(diff,inf) , norm(diff) };

   %    SIRK-2-2
   imethod = imethod + 1;
   mname = 'SIRK-2-2';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  SIRK 2-2' , norm(diff,inf) , norm(diff) };

   %    ESIRK-2-2
   imethod = imethod + 1;
   mname = 'ESIRK-2-2';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' ESIRK 2-2' , norm(diff,inf) , norm(diff) };

   %    Gauss-2-4-IRK
   imethod = imethod + 1;
   mname = 'Gauss-2-4-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' Gauss 2-4' , norm(diff,inf) , norm(diff) };

   %    RadauIIA-2-3-IRK
   imethod = imethod + 1;
   mname = 'RadauIIA-2-3-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  RIIA 2-3' , norm(diff,inf) , norm(diff) };

   %    LobattoIII-2-2-IRK
   imethod = imethod + 1;
   mname = 'LobattoIII-2-2-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  LIII 2-2' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIA-2-2-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIA-2-2-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIA 2-2' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIB-2-2-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIB-2-2-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIB 2-2' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIC-2-2-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIC-2-2-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIC 2-2' , norm(diff,inf) , norm(diff) };

   %    Gauss-3-2-IRK
   imethod = imethod + 1;
   mname = 'Gauss-3-2-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' Gauss 3-2' , norm(diff,inf) , norm(diff) };

   %    RadauI-3-5-IRK
   imethod = imethod + 1;
   mname = 'RadauI-3-5-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '    RI 3-5' , norm(diff,inf) , norm(diff) };

   %    RadauIA-3-5-IRK
   imethod = imethod + 1;
   mname = 'RadauIA-3-5-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '   RIA 3-5' , norm(diff,inf) , norm(diff) };

   %    RadauII-3-5-IRK
   imethod = imethod + 1;
   mname = 'RadauII-3-5-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '   RII 3-5' , norm(diff,inf) , norm(diff) };

   %    RadauIIA-3-5-IRK
   imethod = imethod + 1;
   mname = 'RadauIIA-3-5-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  RIIA 3-5' , norm(diff,inf) , norm(diff) };

   %    LobattoIII-3-4-IRK
   imethod = imethod + 1;
   mname = 'LobattoIII-3-4-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  LIII 3-4' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIA-3-4-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIA-3-4-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIA 3-4' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIB-3-4-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIB-3-4-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIB 3-4' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIC-3-4-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIC-3-4-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIC 3-4' , norm(diff,inf) , norm(diff) };

   %    LobattoIII-4-6-IRK
   imethod = imethod + 1;
   mname = 'LobattoIII-4-6-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  LIII 4-6' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIA-4-6-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIA-4-6-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIA 4-6' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIB-4-6-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIB-4-6-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIB 4-6' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIC-4-6-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIC-4-6-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIC 4-6' , norm(diff,inf) , norm(diff) };

   %    RadauIIA-5-9-IRK
   imethod = imethod + 1;
   mname = 'RadauIIA-5-9-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  RIIA 5-9' , norm(diff,inf) , norm(diff) };

   %    LobattoIII-5-8-IRK
   imethod = imethod + 1;
   mname = 'LobattoIII-5-8-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { '  LIII 5-8' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIA-5-8-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIA-5-8-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIA 5-8' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIB-5-8-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIB-5-8-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIB 5-8' , norm(diff,inf) , norm(diff) };

   %    LobattoIIIC-5-8-IRK
   imethod = imethod + 1;
   mname = 'LobattoIIIC-5-8-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { ' LIIIC 5-8' , norm(diff,inf) , norm(diff) };

   %    Gauss-6-12-IRK
   imethod = imethod + 1;
   mname = 'Gauss-6-12-IRK';
   B = butcher(mname);
   [t,Y,ns] = solve_IRK('f_test', 'J_test', tout, Y0, B, h);
   [ir,ic] = size(Y);  [jr,jc] = size(Ytrue);
   diff = (Y(:,ic)' - Ytrue(jr,:))./Ytrue(jr,:);
   m(imethod,i,1:3) = { 'Gauss 6-12' , norm(diff,inf) , norm(diff) };

end

% set total number of methods to print results on
nmethods = imethod;

% output results for each method
fprintf('\n\n')
for imethod = 1:nmethods
   
   fprintf('  Results for %s:\n',m{imethod,1,1});
   fprintf('      h           errI        err2\n');
   fprintf('  --------------------------------------\n')
   for ih = 1:length(hvals)
      errI(ih) = m{imethod,ih,2};
      err2(ih) = m{imethod,ih,3};
      fprintf('  %10.2e  %10.2e  %10.2e\n',hvals(ih),errI(ih),err2(ih))
   end
   fprintf('  --------------------------------------\n')
   p1 = polyfit(log(hvals),log(errI),1);
   p2 = polyfit(log(hvals),log(err2),1);
   p = (p1(1)+p2(1))/2;
   fprintf('  Measured order of accuracy = %g\n\n\n',p);
   
end


diary off
% end of script
