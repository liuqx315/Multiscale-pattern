function h_test_plots_dirk(probname, pname, ode15s, ...
                           dirk32_c1, dirk32_c2, dirk32_c3, dirk32_c4, dirk32_c5, ...
                           dirk43_c1, dirk43_c2, dirk43_c3, dirk43_c4, dirk43_c5, ...
                           dirk54_c1, dirk54_c2, dirk54_c3, dirk54_c4, dirk54_c5)
% Usage: h_test_plots_dirk(probname, pname, ode15s, ...
%                          dirk32_c1, dirk32_c2, dirk32_c3, dirk32_c4, dirk32_c5, ...
%                          dirk43_c1, dirk43_c2, dirk43_c3, dirk43_c4, dirk43_c5, ...
%                          dirk54_c1, dirk54_c2, dirk54_c3, dirk54_c4, dirk54_c5)
%
%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------

% extract rtol
rtol = ode15s(:,1);

% work vs rtol, DIRK3(2)
figure(1)
loglog(rtol,ode15s(:,3),rtol,dirk32_c2(:,3),rtol,dirk32_c3(:,3),rtol,dirk32_c5(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, DIRK3(2) -- Work vs Tolerance',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG',3)
fname = [ 'work_' , pname , '_dirk32.eps' ];
print('-depsc',fname);

% work vs rtol, DIRK4(3)
figure(2)
loglog(rtol,ode15s(:,3),rtol,dirk43_c2(:,3),rtol,dirk43_c3(:,3),rtol,dirk43_c5(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, DIRK4(3) -- Work vs Tolerance',probname),'FontSize',16)
set(gca,'FontSize',14)
%legend('ode15s','I','PI','PID','iG','eG',3)
legend('ode15s','PI','PID','eG',3)
fname = [ 'work_' , pname , '_dirk43.eps' ];
print('-depsc',fname);

% work vs rtol, DIRK5(4)
figure(3)
loglog(rtol,ode15s(:,3),rtol,dirk54_c2(:,3),rtol,dirk54_c3(:,3),rtol,dirk54_c5(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, DIRK5(4) -- Work vs Tolerance',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG',3)
fname = [ 'work_' , pname , '_dirk54.eps' ];
print('-depsc',fname);



% oversolve, DIRK3(2)
figure(4)
loglog(rtol,ode15s(:,7)./rtol, rtol,dirk32_c2(:,7)./rtol, ...
       rtol,dirk32_c3(:,7)./rtol,  rtol,dirk32_c5(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, DIRK3(2) -- Oversolve',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG',3)
fname = [ 'oversolve_' , pname , '_dirk32.eps' ];
print('-depsc',fname);

% oversolve, DIRK4(3)
figure(5)
loglog(rtol,ode15s(:,7)./rtol, rtol,dirk43_c2(:,7)./rtol, ...
       rtol,dirk43_c3(:,7)./rtol, rtol,dirk43_c5(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, DIRK4(3) -- Oversolve',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG',3)
fname = [ 'oversolve_' , pname , '_dirk43.eps' ];
print('-depsc',fname);

% oversolve, DIRK5(4)
figure(6)
loglog(rtol,ode15s(:,7)./rtol, rtol,dirk54_c2(:,7)./rtol, ...
       rtol,dirk54_c3(:,7)./rtol, rtol,dirk54_c5(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, DIRK5(4) -- Oversolve',probname),'FontSize',16)
set(gca,'FontSize',14)
%legend('ode15s','I','PI','PID','iG','eG',3)
legend('ode15s','PI','PID','eG',3)
fname = [ 'oversolve_' , pname , '_dirk54.eps' ];
print('-depsc',fname);



% work vs err, DIRK3(2)
figure(7)
loglog(ode15s(:,7), ode15s(:,3), dirk32_c2(:,7), dirk32_c2(:,3), ...
       dirk32_c3(:,7), dirk32_c3(:,3), dirk32_c5(:,7), dirk32_c5(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s, DIRK3(2) -- Work vs Error',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG',3)
fname = [ 'work_vs_err_' , pname , '_dirk32.eps' ];
print('-depsc',fname);

% work vs err, DIRK4(3)
figure(8)
loglog(ode15s(:,7), ode15s(:,3), dirk43_c2(:,7), dirk43_c2(:,3), ...
       dirk43_c3(:,7), dirk43_c3(:,3), dirk43_c5(:,7), dirk43_c5(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s, DIRK4(3) -- Work vs Error',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG',3)
fname = [ 'work_vs_err_' , pname , '_dirk43.eps' ];
print('-depsc',fname);

% work vs err, DIRK5(4)
figure(9)
loglog(ode15s(:,7), ode15s(:,3), dirk54_c2(:,7), dirk54_c2(:,3), ...
       dirk54_c3(:,7), dirk54_c3(:,3), dirk54_c5(:,7), dirk54_c5(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s, DIRK5(4) -- Work vs Error',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG',3)
fname = [ 'work_vs_err_' , pname , '_dirk54.eps' ];
print('-depsc',fname);




% work vs rtol, PI
figure(10)
loglog(rtol,ode15s(:,3),rtol,dirk32_c2(:,3),rtol,dirk43_c2(:,3),rtol,dirk54_c2(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Work vs Tolerance (DIRK PI)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'work_' , pname , '_dirk_PI.eps' ];
print('-depsc',fname);

% work vs rtol, PID
figure(11)
loglog(rtol,ode15s(:,3),rtol,dirk32_c3(:,3),rtol,dirk43_c3(:,3),rtol,dirk54_c3(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Work vs Tolerance (DIRK PID)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'work_' , pname , '_dirk_PID.eps' ];
print('-depsc',fname);

% work vs rtol, eG
figure(12)
loglog(rtol,ode15s(:,3),rtol,dirk32_c5(:,3),rtol,dirk43_c5(:,3),rtol,dirk54_c5(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Work vs Tolerance (DIRK eG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'work_' , pname , '_dirk_eG.eps' ];
print('-depsc',fname);



% oversolve, PI
figure(13)
loglog(rtol,ode15s(:,7)./rtol, rtol,dirk32_c2(:,7)./rtol, ...
       rtol,dirk43_c2(:,7)./rtol, rtol,dirk54_c2(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Oversolve (DIRK PI)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'oversolve_' , pname , '_dirk_PI.eps' ];
print('-depsc',fname);

% oversolve, PID
figure(14)
loglog(rtol,ode15s(:,7)./rtol, rtol,dirk32_c3(:,7)./rtol, ...
       rtol,dirk43_c3(:,7)./rtol, rtol,dirk54_c3(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Oversolve (DIRK PID)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'oversolve_' , pname , '_dirk_PID.eps' ];
print('-depsc',fname);

% oversolve, eG
figure(15)
loglog(rtol,ode15s(:,7)./rtol, rtol,dirk32_c5(:,7)./rtol, ...
       rtol,dirk43_c5(:,7)./rtol, rtol,dirk54_c5(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Oversolve (DIRK eG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'oversolve_' , pname , '_dirk_eG.eps' ];
print('-depsc',fname);



% work vs err, PI
figure(16)
loglog(ode15s(:,7), ode15s(:,3), dirk32_c2(:,7), dirk32_c2(:,3), ...
       dirk43_c2(:,7), dirk43_c2(:,3), dirk54_c2(:,7), dirk54_c2(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s -- Work vs Error (DIRK PI)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'work_vs_err_' , pname , '_dirk_PI.eps' ];
print('-depsc',fname);

% work vs err, PID
figure(17)
loglog(ode15s(:,7), ode15s(:,3), dirk32_c3(:,7), dirk32_c3(:,3), ...
       dirk43_c3(:,7), dirk43_c3(:,3), dirk54_c3(:,7), dirk54_c3(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s -- Work vs Error (DIRK PID)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'work_vs_err_' , pname , '_dirk_PID.eps' ];
print('-depsc',fname);

% work vs err, eG
figure(18)
loglog(ode15s(:,7), ode15s(:,3), dirk32_c5(:,7), dirk32_c5(:,3), ...
       dirk43_c5(:,7), dirk43_c5(:,3), dirk54_c5(:,7), dirk54_c5(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s -- Work vs Error (DIRK eG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','DIRK3(2)','DIRK4(3)','DIRK5(4)',3)
fname = [ 'work_vs_err_' , pname , '_dirk_eG.eps' ];
print('-depsc',fname);

