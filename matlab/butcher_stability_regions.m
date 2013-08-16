% Script to plot stability regions for each Butcher table
% comprising ARKode.
%
% Daniel R. Reynolds
% SMU Mathematics
% August 2013
% All Rights Reserved

clear all

% set plotting points
N = 800;
Theta = linspace(0,16*pi,N);

% set list of Butcher tables to plot
tables = {
   'Heun-Euler-ERK',
   'Bogacki-Shampine-ERK',
   'ARK3(2)4L[2]SA-ERK',
   'Zonneveld-4-3-ERK',
   'ARK4(3)6L[2]SA-ERK',
   'Sayfy-Aburub-4-3-ERK',
   'Cash-Karp-ERK',
   'Fehlberg-ERK',
   'Dormand-Prince-ERK',
   'ARK5(4)8L[2]SA-ERK',
   'Verner-6-5-ERK',
   'SDIRK-2-1',
   'Billington-SDIRK',
   'TRBDF2-ESDIRK',
   'Kvaerno(4,2,3)-ESDIRK',
   'ARK3(2)4L[2]SA-ESDIRK',
   'Cash(5,2,4)-SDIRK',
   'Cash(5,3,4)-SDIRK',
   'SDIRK-5-4',
   'Kvaerno(5,3,4)-ESDIRK',
   'ARK4(3)6L[2]SA-ESDIRK',
   'Kvaerno(7,4,5)-ESDIRK',
   'ARK5(4)8L[2]SA-ESDIRK'};

% set list of plot file names
ptype = '-dpng';
pnames = {
   'table0.png',
   'table1.png',
   'table2.png',
   'table3.png',
   'table4.png',
   'table5.png',
   'table6.png',
   'table7.png',
   'table8.png',
   'table9.png',
   'table10.png',
   'table11.png',
   'table12.png',
   'table13.png',
   'table14.png',
   'table15.png',
   'table16.png',
   'table17.png',
   'table18.png',
   'table19.png',
   'table20.png',
   'table21.png',
   'table22.png'};

% iterate over the list of tables
for itable = 1:length(tables)
   % output progress
   table = tables{itable};
   fprintf('Plotting stability region for %s method\n',table);
   % retrieve Butcher table
   B = butcher(table);
   % compute stability region boundary
   [X,Y] = stab_region(B,Theta);
   % create plot
   plot(X,Y), hold on
   xl = xlim();  yl = ylim();
   plot(linspace(xl(1),xl(2),10),zeros(1,10),'b:')
   plot(zeros(1,10),linspace(yl(1),yl(2),10),'b:'), hold off
   xlabel('Re(z)')
   ylabel('Im(z)')
   title(sprintf('Stability boundary for %s method',table))
   print(pnames{itable}, ptype);
end