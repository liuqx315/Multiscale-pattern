function h_test_plots_ark(probname, pname, ode15s, ...
                          ark32_c1, ark32_c2, ark32_c3, ark32_c4, ark32_c5, ark32_c6, ...
                          ark43_c1, ark43_c2, ark43_c3, ark43_c4, ark43_c5, ark43_c6, ...
                          ark54_c1, ark54_c2, ark54_c3, ark54_c4, ark54_c5, ark54_c6)
% Usage: h_test_plots_ark(probname, pname, ode15s, ...
%                          dirk32_c1, dirk32_c2, dirk32_c3, dirk32_c4, dirk32_c5, ...
%                          dirk43_c1, dirk43_c2, dirk43_c3, dirk43_c4, dirk43_c5, ...
%                          dirk54_c1, dirk54_c2, dirk54_c3, dirk54_c4, dirk54_c5)
%
% Daniel R. Reynolds
% SMU Mathematics

% extract rtol
rtol = ode15s(:,1);

% work vs rtol, ARK3(2)
figure(51)
loglog(rtol,ode15s(:,3),rtol,ark32_c2(:,3),rtol,ark32_c3(:,3),...
       rtol,ark32_c5(:,3),rtol,ark32_c6(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, ARK3(2) -- Work vs Tolerance',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'work_' , pname , '_ark32.eps' ];
print('-depsc',fname);

% work vs rtol, ARK4(3)
figure(52)
loglog(rtol,ode15s(:,3),rtol,ark43_c2(:,3),rtol,ark43_c3(:,3),...
       rtol,ark43_c5(:,3),rtol,ark43_c6(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, ARK4(3) -- Work vs Tolerance',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'work_' , pname , '_ark43.eps' ];
print('-depsc',fname);

% work vs rtol, ARK5(4)
figure(53)
loglog(rtol,ode15s(:,3),rtol,ark54_c2(:,3),rtol,ark54_c3(:,3),...
       rtol,ark54_c5(:,3),rtol,ark54_c6(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, ARK5(4) -- Work vs Tolerance',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'work_' , pname , '_ark54.eps' ];
print('-depsc',fname);



% oversolve, ARK3(2)
figure(54)
loglog(rtol,ode15s(:,7)./rtol, rtol,ark32_c2(:,7)./rtol, ...
       rtol,ark32_c3(:,7)./rtol, rtol,ark32_c5(:,7)./rtol, ...
       rtol,ark32_c6(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, ARK3(2) -- Oversolve',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'oversolve_' , pname , '_ark32.eps' ];
print('-depsc',fname);

% oversolve, ARK4(3)
figure(55)
loglog(rtol,ode15s(:,7)./rtol, rtol,ark43_c2(:,7)./rtol, ...
       rtol,ark43_c3(:,7)./rtol, ...
       rtol,ark43_c5(:,7)./rtol, rtol,ark43_c6(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, ARK4(3) -- Oversolve',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'oversolve_' , pname , '_ark43.eps' ];
print('-depsc',fname);

% oversolve, ARK5(4)
figure(56)
loglog(rtol,ode15s(:,7)./rtol, rtol,ark54_c2(:,7)./rtol, ...
       rtol,ark54_c3(:,7)./rtol, ...
       rtol,ark54_c5(:,7)./rtol, rtol,ark54_c6(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s, ARK5(4) -- Oversolve',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'oversolve_' , pname , '_ark54.eps' ];
print('-depsc',fname);



% work vs err, ARK3(2)
figure(57)
loglog(ode15s(:,7), ode15s(:,3), ark32_c2(:,7), ark32_c2(:,3), ...
       ark32_c3(:,7), ark32_c3(:,3), ...
       ark32_c5(:,7), ark32_c5(:,3) ,ark32_c6(:,7), ark32_c6(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s, ARK3(2) -- Work vs Error',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'work_vs_err_' , pname , '_ark32.eps' ];
print('-depsc',fname);

% work vs err, ARK4(3)
figure(58)
loglog(ode15s(:,7), ode15s(:,3), ark43_c2(:,7), ark43_c2(:,3), ...
       ark43_c3(:,7), ark43_c3(:,3), ...
       ark43_c5(:,7), ark43_c5(:,3) ,ark43_c6(:,7), ark43_c6(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s, ARK4(3) -- Work vs Error',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'work_vs_err_' , pname , '_ark43.eps' ];
print('-depsc',fname);

% work vs err, ARK5(4)
figure(59)
loglog(ode15s(:,7), ode15s(:,3), ark54_c2(:,7), ark54_c2(:,3), ...
       ark54_c3(:,7), ark54_c3(:,3), ...
       ark54_c5(:,7), ark54_c5(:,3) ,ark54_c6(:,7), ark54_c6(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s, ARK5(4) -- Work vs Error',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','PI','PID','eG','ieG',3)
fname = [ 'work_vs_err_' , pname , '_ark54.eps' ];
print('-depsc',fname);



% work vs rtol, PI
figure(60)
loglog(rtol,ode15s(:,3),rtol,ark32_c2(:,3),rtol,ark43_c2(:,3),...
       rtol,ark54_c2(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Work vs Tolerance (ARK PI)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'work_' , pname , '_ark_PI.eps' ];
print('-depsc',fname);

figure(61)
loglog(rtol,ode15s(:,3),rtol,ark32_c3(:,3),rtol,ark43_c3(:,3),...
       rtol,ark54_c3(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Work vs Tolerance (ARK PID)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'work_' , pname , '_ark_PID.eps' ];
print('-depsc',fname);

figure(62)
loglog(rtol,ode15s(:,3),rtol,ark32_c5(:,3),rtol,ark43_c5(:,3),...
       rtol,ark54_c5(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Work vs Tolerance (ARK eG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'work_' , pname , '_ark_eG.eps' ];
print('-depsc',fname);

figure(63)
loglog(rtol,ode15s(:,3),rtol,ark32_c6(:,3),rtol,ark43_c6(:,3),...
       rtol,ark54_c6(:,3),'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Work vs Tolerance (ARK ieG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'work_' , pname , '_ark_ieG.eps' ];
print('-depsc',fname);



% oversolve, PI
figure(64)
loglog(rtol,ode15s(:,7)./rtol, rtol,ark32_c2(:,7)./rtol, ...
       rtol,ark43_c2(:,7)./rtol, rtol,ark54_c2(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Oversolve (ARK PI)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'oversolve_' , pname , '_ark_PI.eps' ];
print('-depsc',fname);

figure(65)
loglog(rtol,ode15s(:,7)./rtol, rtol,ark32_c3(:,7)./rtol, ...
       rtol,ark43_c3(:,7)./rtol, rtol,ark54_c3(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Oversolve (ARK PID)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'oversolve_' , pname , '_ark_PID.eps' ];
print('-depsc',fname);

figure(66)
loglog(rtol,ode15s(:,7)./rtol, rtol,ark32_c5(:,7)./rtol, ...
       rtol,ark43_c5(:,7)./rtol, rtol,ark54_c5(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Oversolve (ARK eG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'oversolve_' , pname , '_ark_eG.eps' ];
print('-depsc',fname);

figure(67)
loglog(rtol,ode15s(:,7)./rtol, rtol,ark32_c6(:,7)./rtol, ...
       rtol,ark43_c6(:,7)./rtol, rtol,ark54_c6(:,7)./rtol, 'LineWidth',2)
xlabel('r_{tol}','FontSize',14)
title(sprintf('%s -- Oversolve (ARK ieG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'oversolve_' , pname , '_ark_ieG.eps' ];
print('-depsc',fname);



% work vs err, PI
figure(68)
loglog(ode15s(:,7), ode15s(:,3), ark32_c2(:,7), ark32_c2(:,3), ...
       ark43_c2(:,7), ark43_c2(:,3), ark54_c2(:,7), ark54_c2(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s -- Work vs Error (ARK PI)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'work_vs_err_' , pname , '_ark_PI.eps' ];
print('-depsc',fname);

figure(69)
loglog(ode15s(:,7), ode15s(:,3), ark32_c3(:,7), ark32_c3(:,3), ...
       ark43_c3(:,7), ark43_c3(:,3), ark54_c3(:,7), ark54_c3(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s -- Work vs Error (ARK PID)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'work_vs_err_' , pname , '_ark_PID.eps' ];
print('-depsc',fname);

figure(70)
loglog(ode15s(:,7), ode15s(:,3), ark32_c5(:,7), ark32_c5(:,3), ...
       ark43_c5(:,7), ark43_c5(:,3), ark54_c5(:,7), ark54_c5(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s -- Work vs Error (ARK eG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'work_vs_err_' , pname , '_ark_eG.eps' ];
print('-depsc',fname);

figure(71)
loglog(ode15s(:,7), ode15s(:,3), ark32_c6(:,7), ark32_c6(:,3), ...
       ark43_c6(:,7), ark43_c6(:,3), ark54_c6(:,7), ark54_c6(:,3), 'LineWidth',2)
xlabel('error','FontSize',14)
ylabel('work','FontSize',14)
title(sprintf('%s -- Work vs Error (ARK ieG)',probname),'FontSize',16)
set(gca,'FontSize',14)
legend('ode15s','ARK3(2)','ARK4(3)','ARK5(4)',3)
fname = [ 'work_vs_err_' , pname , '_ark_ieG.eps' ];
print('-depsc',fname);

