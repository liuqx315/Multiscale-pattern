%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2013, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------
% Script to plot stability regions for each Butcher table
% comprising ARKode.

clear all

% set plotting points
N = 10000;
Theta = linspace(0,16*pi,N);

% set graphics type
gtype = 'png';


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
ptype = [ '-d' , gtype ];
pnames = {
   'stab_region_0',
   'stab_region_1',
   'stab_region_2',
   'stab_region_3',
   'stab_region_4',
   'stab_region_5',
   'stab_region_6',
   'stab_region_7',
   'stab_region_8',
   'stab_region_9',
   'stab_region_10',
   'stab_region_11',
   'stab_region_12',
   'stab_region_13',
   'stab_region_14',
   'stab_region_15',
   'stab_region_16',
   'stab_region_17',
   'stab_region_18',
   'stab_region_19',
   'stab_region_20',
   'stab_region_21',
   'stab_region_22'};

% iterate over the list of tables
for itable = 1:length(tables)
   % output progress
   table = tables{itable};
   fprintf('Plotting stability region for %s method\n',table);
   % retrieve Butcher table and extract components
   B = butcher(table);
   [m,n] = size(B);
   s = n-1;
   b = B(s+1,2:n)';
   A = B(1:s,2:n);
   b2 = B(s+2,2:n)';
   % compute method stability region boundary and create plot
   [X,Y] = stab_region(A,b,Theta);
   plot(X,Y,'b-'), hold on
   % compute embedding stability region boundary, add to plot
   [X,Y] = stab_region(A,b2,Theta);
   plot(X,Y,'r-')
   xl = xlim();  yl = ylim();
   plot(linspace(xl(1),xl(2),10),zeros(1,10),'b:')
   plot(zeros(1,10),linspace(yl(1),yl(2),10),'b:'), hold off
   xlabel('Re(z)')
   ylabel('Im(z)')
   title(sprintf('Stability boundary for %s method',table))
   legend('method','embedding')
   pname = [ pnames{itable}, '.', gtype ];
   print(pname, ptype);
end