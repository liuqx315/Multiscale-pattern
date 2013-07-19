:tocdepth: 3


.. _Mathematics.Butcher:

Butcher tables
======================

Here we catalogue the full set of Butcher tables included in ARKode.
We group these into three categories: *explicit*, *implicit* and
*additive*.  However, since the methods that comprise an additive
Runge Kutta method are themselves explicit and implicit, their
component Butcher tables are listed both within their separate
sections, but are referenced together in the additive section.

In each of the following tables, we use the following notation:

.. math::

   \begin{array}{r|ccc}
     c_1 & a_{1,1} & a_{1,2} & a_{1,3} \\
     c_2 & a_{2,1} & a_{2,2} & a_{2,3} \\
     c_3 & a_{3,1} & a_{3,2} & a_{3,3} \\
     \hline
     q & b_1 & b_2 & b_3 \\
     p & \tilde{b}_1 & \tilde{b}_2 & \tilde{b}_3
   \end{array}

where here the method and embedding share stage :math:`A` and
:math:`c` values, but use their stages :math:`z_i` differently through
the coefficients :math:`b` and :math:`\tilde{b}` to generate methods
of orders :math:`q` (the main method) and :math:`p` (the embedding,
typically :math:`q = p+1`).


.. _Mathematics.Butcher.explicit:

Explicit Butcher tables
---------------------------

In the category of explicit Runge-Kutta methods, ARKode includes
methods that have orders 2 through 6, with embeddings that are of
orders 1 through 5.


.. _Mathematics.Butcher.Heun_Euler:

Heun-Euler
^^^^^^^^^^^^

Butcher table number 0 for :func:`ARKodeSetERKTableNum()`.  This is
the default 2nd order explicit method.

.. math::

   \begin{array}{r|cc}
     0 & 0 & 0 \\
     1 & 1 & 0 \\
     \hline
     2 & 1/2 & 1/2 \\
     1 & 1 & 0
   \end{array}



.. _Mathematics.Butcher.Bogacki_Shampine:

Bogacki-Shampine
^^^^^^^^^^^^^^^^^^

Butcher table number 1 for :func:`ARKodeSetERKTableNum()`.  This is
the default 3rd order explicit method.

.. math::

   \begin{array}{r|cccc}
     0 &   0 & 0 & 0 & 0 \\
     1/2 & 1/2 & 0 & 0 & 0 \\
     3/4 & 0 & 3/4 & 0 & 0 \\
     1   & 2/9 & 1/3 & 4/9 & 0 \\
     \hline
     3 & 2/9 & 1/3 & 4/9 \\
     2 & 7/24 & 1/4 & 1/3 & 1/8
   \end{array}


.. _Mathematics.Butcher.ARK_3_2_E:

ARK3(2)4L[2]SA (explicit)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 2 for :func:`ARKodeSetERKTableNum()`.

.. math::

   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     \frac{1767732205903}{2027836641118} & \frac{1767732205903}{2027836641118} & 0 & 0 & 0 \\
     3/5 & \frac{5535828885825}{10492691773637} & \frac{788022342437}{10882634858940} & 0 & 0 \\
     1 & \frac{6485989280629}{16251701735622} & -\frac{4246266847089}{9704473918619} & \frac{10755448449292}{10357097424841} & 0 \\
     \hline
     3 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     2 & \frac{2756255671327}{12835298489170} & -\frac{10771552573575}{22201958757719} & \frac{9247589265047}{10645013368117} & \frac{2193209047091}{5459859503100}
   \end{array}


.. _Mathematics.Butcher.Zonneveld:

Zonneveld
^^^^^^^^^^^^

Butcher table number 3 for :func:`ARKodeSetERKTableNum()`.  This is
the default 4th order explicit method.

.. math::

   \begin{array}{r|ccccc}
       0 & 0 & 0 & 0 & 0 & 0 \\
     1/2 & 1/2 & 0 & 0 & 0 & 0 \\
     1/2 & 0 & 1/2 & 0 & 0 & 0 \\
       1 & 0 & 0 & 1 & 0 & 0 \\
     3/4 & 5/32 & 7/32 & 13/32 & -1/32 & 0 \\
     \hline
     4 & 1/6 & 1/3 & 1/3 & 1/6 & 0 \\
     3 & -1/2 & 7/3 & 7/3 & 13/6 & -16/3
   \end{array}


.. _Mathematics.Butcher.ARK_4_3_E:

ARK4(3)6L[2]SA (explicit)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 4 for :func:`ARKodeSetERKTableNum()`.

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac12 & \frac12 & 0 & 0 & 0 & 0 & 0 \\
     \frac{83}{250} & \frac{13861}{62500} & \frac{6889}{62500} & 0 & 0 & 0 & 0 \\
     \frac{31}{50} & -\frac{116923316275}{2393684061468} & -\frac{2731218467317}{15368042101831} & \frac{9408046702089}{11113171139209} & 0 & 0 & 0 \\
     \frac{17}{20} & -\frac{451086348788}{2902428689909} & -\frac{2682348792572}{7519795681897} & \frac{12662868775082}{11960479115383} & \frac{3355817975965}{11060851509271} & 0 & 0 \\
     1 & \frac{647845179188}{3216320057751} & \frac{73281519250}{8382639484533} & \frac{552539513391}{3454668386233} & \frac{3354512671639}{8306763924573} & \frac{4040}{17871} & 0 \\
     \hline
     4 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac14 \\
     3 & \frac{4586570599}{29645900160} & 0 & \frac{178811875}{945068544} & \frac{814220225}{1159782912} & -\frac{3700637}{11593932} & \frac{61727}{225920}
   \end{array}


.. _Mathematics.Butcher.Sayfy_Aburub:

Sayfy-Aburub-4-3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 5 for :func:`ARKodeSetERKTableNum()`.

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     1/2 & 1/2 & 0 & 0 & 0 & 0 & 0 \\
     1 & -1 & 2 & 0 & 0 & 0 & 0 \\
     1 & 1/6 & 2/3 & 1/6 & 0 & 0 & 0 \\
     1/2 & 0.137 & 0.226 & 0.173 & 0 & 0 & 0 \\
     1 & 0.452 & -0.904 & -0.548 & 0 & 2 & 0 \\
     \hline
     4 & 1/6 & 1/3 & 1/12 & 0 & 1/3 & 1/12 \\
     3 & 1/6 & 2/3 & 1/6 & 0 & 0 & 0
   \end{array}



.. _Mathematics.Butcher.Cash-Karp:

Cash-Karp
^^^^^^^^^^^^

Butcher table number 6 for :func:`ARKodeSetERKTableNum()`.  This is
the default 5th order explicit method.

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     1/5 & 1/5 & 0 & 0 & 0 & 0 & 0 \\
     3/10 & 3/40 & 9/40 & 0 & 0 & 0 & 0 \\
     3/5 & 3/10 & -9/10 & 6/5 & 0 & 0 & 0 \\
     1 & -11/54 & 5/2 & -70/27 & 35/27 & 0 & 0 \\
     7/8 & 1631/55296 & 175/512 & 575/13824 & 44275/110592 & 253/4096 & 0 \\
     \hline
     5 & 37/348 & 0 & 250/621 & 125/594 & 0 & 512/1771 \\
     4 & 2825/27648 & 0 & 18575/48384 & 13525/55296 & 277/14336 & 1/4
   \end{array}




.. _Mathematics.Butcher.Fehlberg:

Fehlberg
^^^^^^^^^^^^

Butcher table number 7 for :func:`ARKodeSetERKTableNum()`.

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     1/4 & 1/4 & 0 & 0 & 0 & 0 & 0 \\
     3/8 & 3/32 & 9/32 & 0 & 0 & 0 & 0 \\
     12/13 & 1932/2197 & -7200/2197 & 7296/2197 & 0 & 0 & 0 \\
     1 & 439/216 & -8 & 3680/513 & -845/4104 & 0 & 0 \\
     1/2 & -8/27 & 2 & -3544/2565 & 1859/4104 & -11/40 & 0 \\
     \hline
     5 & 25/216 & 0 & 1408/2565 & 2197/4104 & -1/5 & 0 \\
     4 & 16/135 & 0 & 6656/12825 & 28561/56430 & -9/50 & 2/55
   \end{array}



.. _Mathematics.Butcher.Dormand_Prince:

Dormand-Prince
^^^^^^^^^^^^^^^^^

Butcher table number 8 for :func:`ARKodeSetERKTableNum()`.

.. math::

   \begin{array}{r|ccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     1/5 & 1/5 & 0 & 0 & 0 & 0 & 0 & 0 \\
     3/10 & 3/40 & 9/40 & 0 & 0 & 0 & 0 & 0 \\
     4/5 & 44/45 & -56/15 & 32/9 & 0 & 0 & 0 & 0 \\
     8/9 & 19372/6561 & -25360/2187 & 64448/6561 & -212/729 & 0 & 0 & 0 \\
     1 & 9017/3168 & -355/33 & 46732/5247 & 49/176 & -5103/18656 & 0 & 0 \\
     1 & 35/384 & 0 & 500/1113 & 125/192 & -2187/6784 & 11/84 & 0 \\
     \hline
     5 & 35/384 & 0 & 500/1113 & 125/192 & -2187/6784 & 11/84 & 0 \\
     4 & 5179/57600 & 0 & 7571/16695 & 393/640 & -92097/339200 & 187/2100 & 1/40
   \end{array}



.. _Mathematics.Butcher.ARK_5_4_E:

ARK5(4)8L[2]SA (explicit)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 9 for :func:`ARKodeSetERKTableNum()`.

.. math::

   \begin{array}{r|cccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{41}{100} & \frac{41}{100} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{2935347310677}{11292855782101} & \frac{367902744464}{2072280473677} & \frac{677623207551}{8224143866563} & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1426016391358}{7196633302097} & \frac{1268023523408}{10340822734521} & 0 & \frac{1029933939417}{13636558850479} & 0 & 0 & 0 & 0 & 0 \\
     \frac{92}{100} & \frac{14463281900351}{6315353703477} & 0 & \frac{66114435211212}{5879490589093} & -\frac{54053170152839}{4284798021562} & 0 & 0 & 0 & 0 \\
     \frac{24}{100} & \frac{14090043504691}{34967701212078} & 0 & \frac{15191511035443}{11219624916014} & -\frac{18461159152457}{12425892160975} & -\frac{281667163811}{9011619295870} & 0 & 0 & 0 \\
     \frac{3}{5} & \frac{19230459214898}{13134317526959} & 0 & \frac{21275331358303}{2942455364971} & -\frac{38145345988419}{4862620318723} & -\frac{1}{8} & -\frac{1}{8} & 0 & 0 \\
     1 & -\frac{19977161125411}{11928030595625} & 0 & -\frac{40795976796054}{6384907823539} & \frac{177454434618887}{12078138498510} & \frac{782672205425}{8267701900261} & -\frac{69563011059811}{9646580694205} & \frac{7356628210526}{4942186776405} & 0 \\
     \hline
     5 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
     4 & -\frac{975461918565}{9796059967033} & 0 & 0 & \frac{78070527104295}{32432590147079} & -\frac{548382580838}{3424219808633} & -\frac{33438840321285}{15594753105479} & \frac{3629800801594}{4656183773603} & \frac{4035322873751}{18575991585200}
   \end{array}




.. _Mathematics.Butcher.Verner-6-5:

Verner-6-5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 10 for :func:`ARKodeSetERKTableNum()`.  This is
the default 6th order explicit method.

.. math::

   \begin{array}{r|cccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     1/6 & 1/6 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     4/15 & 4/75 & 16/75 & 0 & 0 & 0 & 0 & 0 & 0 \\
     2/3 & 5/6 & -8/3 & 5/2 & 0 & 0 & 0 & 0 & 0 \\
     5/6 & -165/64 & 55/6 & -425/64 & 85/96 & 0 & 0 & 0 & 0 \\
     1 & 12/5 & -8 & 4015/612 & -11/36 & 88/255 & 0 & 0 & 0 \\
     1/15 & -8263/15000 & 124/75 & -643/680 & -81/250 & 2484/10625 & 0 & 0 & 0 \\
     1 & 3501/1720 & -300/43 & 297275/52632 & -319/2322 & 24068/84065 & 0 & 3850/26703 & 0 \\
     \hline
     6 & 3/40 & 0 & 875/2244 & 23/72 & 264/1955 & 0 & 125/11592 & 43/616 \\
     5 & 13/160 & 0 & 2375/5984 & 5/16 & 12/85 & 3/44 & 0 & 0
   \end{array}




.. _Mathematics.Butcher.implicit:

Implicit Butcher tables
---------------------------


In the category of diagonally implicit Runge-Kutta methods, ARKode
includes methods that have orders 2 through 5, with embeddings that are of
orders 1 through 4.


.. _Mathematics.Butcher.SDIRK-2-1:

SDIRK 2(1)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 11 for :func:`ARKodeSetIRKTableNum()`.  This is
the default 2nd order implicit method.

.. math::

   \begin{array}{r|cc}
     0 & 1 & 0 \\
     1 & -1 & 1 \\
     \hline
     2 & 1/2 & 1/2 \\
     1 & 1 & 0
   \end{array}



.. _Mathematics.Butcher.Billington:

Billington
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 12 for :func:`ARKodeSetIRKTableNum()`.

.. math::

   \begin{array}{r|ccc}
     0.292893218813 & 0.292893218813 & 0 & 0 \\
     1.091883092037 & 0.798989873223 & 0.292893218813 & 0 \\
     1.292893218813 & 0.740789228841 & 0.259210771159 & 0.292893218813 \\
     \hline
     3 & 0.691665115992 & 0.503597029883 & -0.195262145876 \\
     2 & 0.740789228840 & 0.259210771159 & 0
   \end{array}




.. _Mathematics.Butcher.TRBDF2:

TRBDF2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 13 for :func:`ARKodeSetIRKTableNum()`.

.. math::

   \begin{array}{r|ccc}
     0 & 0 & 0 & 0 \\
     2-\sqrt{2} & \frac{2-\sqrt{2}}{2} & \frac{2-\sqrt{2}}{2} & 0 \\
     1 & \frac{\sqrt{2}}{4} & \frac{\sqrt{2}}{4} & \frac{2-\sqrt{2}}{2} \\
     \hline
     3 & \frac{1-\frac{\sqrt{2}}{4}}{3} & \frac{\frac{3\sqrt{2}}{4}+1}{3} & \frac{2-\sqrt{2}}{6} \\
     2 & \frac{\sqrt{2}}{4} & \frac{\sqrt{2}}{4} & \frac{2-\sqrt{2}}{2}
   \end{array}




.. _Mathematics.Butcher.Kvaerno_4_2_3:

Kvaerno(4,2,3)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 14 for :func:`ARKodeSetIRKTableNum()`.  This is
the default 3rd order implicit method.

.. math::

   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     0.871733043 & 0.4358665215 & 0.4358665215 & 0 & 0 \\
     1 & 0.490563388419108 & 0.073570090080892 & 0.4358665215 & 0 \\
     1 & 0.308809969973036 & 1.490563388254106 & -1.235239879727145 & 0.4358665215 \\
     \hline
     3 & 0.308809969973036 & 1.490563388254106 & -1.235239879727145 & 0.4358665215 \\
     2 & 0.490563388419108 & 0.073570090080892 & 0.4358665215 & 0
   \end{array}



.. _Mathematics.Butcher.ARK_3_2_I:

ARK3(2)4L[2]SA (implicit)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 15 for :func:`ARKodeSetIRKTableNum()`.

.. math::

   \begin{array}{r|cccc}
     0 & 0 & 0 & 0 & 0 \\
     \frac{1767732205903}{2027836641118} & \frac{1767732205903}{4055673282236} & \frac{1767732205903}{4055673282236} & 0 & 0 \\
     \frac{3}{5} & \frac{2746238789719}{10658868560708} & -\frac{640167445237}{6845629431997} & \frac{1767732205903}{4055673282236} & 0 \\
     1 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     \hline
     3 & \frac{1471266399579}{7840856788654} & -\frac{4482444167858}{7529755066697} & \frac{11266239266428}{11593286722821} & \frac{1767732205903}{4055673282236} \\
     2 & \frac{2756255671327}{12835298489170} & -\frac{10771552573575}{22201958757719} & \frac{9247589265047}{10645013368117} & \frac{2193209047091}{5459859503100}
   \end{array}




.. _Mathematics.Butcher.Cash_5_2_4:

Cash(5,2,4)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 16 for :func:`ARKodeSetIRKTableNum()`.

.. math::

   \begin{array}{r|ccccc}
     0.435866521508 & 0.435866521508 & 0 & 0 & 0 & 0 \\
     -0.7 & -1.13586652150 & 0.435866521508 & 0 & 0 & 0 \\
     0.8 & 1.08543330679 & -0.721299828287 & 0.435866521508 & 0 & 0 \\
     0.924556761814 & 0.416349501547 & 0.190984004184 & -0.118643265417 & 0.435866521508 & 0 \\
     1 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     \hline
     4 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     2 & 1.05646216107052 & -0.0564621610705236 & 0 & 0 & 0
   \end{array}




.. _Mathematics.Butcher.Cash_5_3_4:

Cash(5,3,4)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 17 for :func:`ARKodeSetIRKTableNum()`.

.. math::

   \begin{array}{r|ccccc}
     0.435866521508 & 0.435866521508 & 0 & 0 & 0 & 0 \\
     -0.7 & -1.13586652150 & 0.435866521508 & 0 & 0 & 0 \\
     0.8 & 1.08543330679 & -0.721299828287 & 0.435866521508 & 0 & 0 \\
     0.924556761814 & 0.416349501547 & 0.190984004184 & -0.118643265417 & 0.435866521508 & 0 \\
     1 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     \hline
     4 & 0.896869652944 & 0.0182725272734 & -0.0845900310706 & -0.266418670647 & 0.435866521508 \\
     3 & 0.776691932910 & 0.0297472791484 & -0.0267440239074 & 0.220304811849 & 0
   \end{array}



.. _Mathematics.Butcher.SDIRK-5-4:

SDIRK 5(4)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 18 for :func:`ARKodeSetIRKTableNum()`.  This is
the default 4th order implicit method.

.. math::

   \begin{array}{r|ccccc}
     1/4 & 1/4 & 0 & 0 & 0 & 0 \\
     3/4 & 1/2 & 1/4 & 0 & 0 & 0 \\
     11/20 & 17/50 & -1/25 & 1/4 & 0 & 0 \\
     1/2 & 371/1360 & -137/2720 & 15/544 & 1/4 & 0 \\
     1 & 25/24 & -49/48 & 125/16 & -85/12 & 1/4 \\
     \hline
     4 & 25/24 & -49/48 & 125/16 & -85/12 & 1/4 \\
     3 & 59/48 & -17/96 & 225/32 & -85/12 & 0
   \end{array}





.. _Mathematics.Butcher.Kvaerno_5_3_4:

Kvaerno(5,3,4)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 19 for :func:`ARKodeSetIRKTableNum()`.

.. math::

   \begin{array}{r|ccccc}
     0 & 0 & 0 & 0 & 0 & 0 \\
     0.871733043 & 0.4358665215  & 0.4358665215  & 0 & 0 & 0 \\
     0.468238744853136 & 0.140737774731968 & -0.108365551378832 & 0.4358665215 & 0 & 0 \\
     1 & 0.102399400616089 & -0.376878452267324 & 0.838612530151233 & 0.4358665215 & 0 \\
     1 & 0.157024897860995 & 0.117330441357768 & 0.61667803039168 & -0.326899891110444 & 0.4358665215 \\
     \hline
     4 & 0.157024897860995 & 0.117330441357768 & 0.61667803039168 & -0.326899891110444 & 0.4358665215 \\
     3 & 0.102399400616089 & -0.376878452267324 & 0.838612530151233 & 0.4358665215 & 0
   \end{array}




.. _Mathematics.Butcher.ARK_4_3_I:

ARK4(3)6L[2]SA (implicit)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 20 for :func:`ARKodeSetIRKTableNum()`.

.. math::

   \begin{array}{r|cccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{1}{2} & \frac{1}{4} & \frac{1}{4} & 0 & 0 & 0 & 0 \\
     \frac{83}{250} & \frac{8611}{62500} & -\frac{1743}{31250} & \frac{1}{4} & 0 & 0 & 0 \\
     \frac{31}{50} & \frac{5012029}{34652500} & -\frac{654441}{2922500} & \frac{174375}{388108} & \frac{1}{4} & 0 & 0 \\
     \frac{17}{20} & \frac{15267082809}{155376265600} & -\frac{71443401}{120774400} & \frac{730878875}{902184768} & \frac{2285395}{8070912} & \frac{1}{4} & 0 \\
     1 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac{1}{4} \\
     \hline
     4 & \frac{82889}{524892} & 0 & \frac{15625}{83664} & \frac{69875}{102672} & -\frac{2260}{8211} & \frac{1}{4} \\
     3 & \frac{4586570599}{29645900160} & 0 & \frac{178811875}{945068544} & \frac{814220225}{1159782912} & -\frac{3700637}{11593932} & \frac{61727}{225920}
   \end{array}




.. _Mathematics.Butcher.Kvaerno_7_4_5:

Kvaerno(7,4,5)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 21 for :func:`ARKodeSetIRKTableNum()`.

.. math::

   \begin{array}{r|ccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     0.52 & 0.26 & 0.26 & 0 & 0 & 0 & 0 & 0 \\
     1.230333209967908 & 0.13 & 0.84033320996790809 & 0.26 & 0 & 0 & 0 & 0 \\
     0.895765984350076 & 0.22371961478320505 & 0.47675532319799699 & -0.06470895363112615 & 0.26 & 0 & 0 & 0 \\
     0.436393609858648 & 0.16648564323248321 & 0.10450018841591720 & 0.03631482272098715 & -0.13090704451073998 & 0.26 & 0 & 0 \\
     1 & 0.13855640231268224 & 0 & -0.04245337201752043 & 0.02446657898003141 & 0.61943039072480676 & 0.26 & 0 \\
     1 & 0.13659751177640291 & 0 & -0.05496908796538376 & -0.04118626728321046 & 0.62993304899016403 & 0.06962479448202728 & 0.26 \\
     \hline
     5 & 0.13659751177640291 & 0 & -0.05496908796538376 & -0.04118626728321046 & 0.62993304899016403 & 0.06962479448202728 & 0.26 \\
     4 & 0.13855640231268224 & 0 & -0.04245337201752043 & 0.02446657898003141 & 0.61943039072480676 & 0.26 & 0
   \end{array}




.. _Mathematics.Butcher.ARK_5_4_I:

ARK5(4)8L[2]SA (implicit)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Butcher table number 22 for :func:`ARKodeSetIRKTableNum()`.  This is
the default 5th order implicit method.

.. math::

   \begin{array}{r|cccccccc}
     0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{41}{100} & \frac{41}{200} & \frac{41}{200} & 0 & 0 & 0 & 0 & 0 & 0 \\
     \frac{2935347310677}{11292855782101} & \frac{41}{400} & -\frac{567603406766}{11931857230679} & \frac{41}{200} & 0 & 0 & 0 & 0 & 0 \\
     \frac{1426016391358}{7196633302097} & \frac{683785636431}{9252920307686} & 0 & -\frac{110385047103}{1367015193373} & \frac{41}{200} & 0 & 0 & 0 & 0 \\
     \frac{92}{100} & \frac{3016520224154}{10081342136671} & 0 & \frac{30586259806659}{12414158314087} & -\frac{22760509404356}{11113319521817} & \frac{41}{200} & 0 & 0 & 0 \\
     \frac{24}{100} & \frac{218866479029}{1489978393911} & 0 & \frac{638256894668}{5436446318841} & -\frac{1179710474555}{5321154724896} & -\frac{60928119172}{8023461067671} & \frac{41}{200} & 0 & 0 \\
     \frac{3}{5} & \frac{1020004230633}{5715676835656} & 0 & \frac{25762820946817}{25263940353407} & -\frac{2161375909145}{9755907335909} & -\frac{211217309593}{5846859502534} & -\frac{4269925059573}{7827059040749} & \frac{41}{200} & 0 \\
     1 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
     \hline
     5 & -\frac{872700587467}{9133579230613} & 0 & 0 & \frac{22348218063261}{9555858737531} & -\frac{1143369518992}{8141816002931} & -\frac{39379526789629}{19018526304540} & \frac{32727382324388}{42900044865799} & \frac{41}{200} \\
     4 & -\frac{975461918565}{9796059967033} & 0 & 0 & \frac{78070527104295}{32432590147079} & -\frac{548382580838}{3424219808633} & -\frac{33438840321285}{15594753105479} & \frac{3629800801594}{4656183773603} & \frac{4035322873751}{18575991585200}
   \end{array}






.. _Mathematics.Butcher.additive:

Additive Butcher tables
---------------------------

In the category of additive Runge-Kutta methods for split implicit and
explicit calculations, ARKode includes methods that have orders 3
through 5, with embeddings that are of orders 2 through 4.  These
Butcher table pairs are as follows:

* 3rd-order pair:  :ref:`Mathematics.Butcher.ARK_3_2_E` with
  :ref:`Mathematics.Butcher.ARK_3_2_I`, corresponding to Butcher
  tables 2 and 15 for :func:`ARKodeSetARKTableNum()`.

* 4th-order pair:  :ref:`Mathematics.Butcher.ARK_4_3_E` with
  :ref:`Mathematics.Butcher.ARK_4_3_I`, corresponding to Butcher 
  tables 4 and 20 for :func:`ARKodeSetARKTableNum()`. 

* 5th-order pair:  :ref:`Mathematics.Butcher.ARK_5_4_E` with
  :ref:`Mathematics.Butcher.ARK_5_4_I`, corresponding to Butcher
  tables 9 and 22 for :func:`ARKodeSetARKTableNum()`

