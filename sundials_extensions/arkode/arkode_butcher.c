/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 This is the implementation file for the known Butcher tables 
 for the ARKODE solver.

 VERIFIED TO CORRECTLY FILL IN TABLES!
---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


#define SQRT2 RSqrt(RCONST(2.0))

/*---------------------------------------------------------------
 Returns butcher table for Runge Kutta methods.  
 Input:
     imeth -- integer key for the desired method (see below)
 Outputs:
     s -- int number of stages for the requested method [s <= ARK_S_MAX = 8]
     q -- integer, theoretical order of accuracy for the method
     p -- integer, theoretical order of accuracy for the embedding
     A[s][s] -- realtype Butcher table coefficients
     b[s] -- realtype root node coefficients
     c[s] -- realtype canopy node coefficients
     b2[s] -- realtype embedding coefficients
     bdense[s][s] -- realtype dense output coefficients, only 
               filled in if known.

 Allowed 'method' names and properties:

   imeth    name             type   s    q    p  dense  A  L 
  -----------------------------------------------------------
    0     Heun-Euler          ERK   2    2    1   no    -  - 
    1     Billington        SDIRK   3    3    2   no    -  - 
    2     TRBDF2           ESDIRK   3    3    2   no    -  - 
    3     TRX2             ESDIRK   3    4    2   no    -  - 
    4     Bogacki-Shampine    ERK   4    4?   3?  no    -  - 
    5     ARK3(2)4L[2]SA      ERK   4    3    2   yes   -  - 
    6     ARK3(2)4L[2]SA   ESDIRK   4    3    2   yes   X  X 
    7     SDIRK-5-4         SDIRK   5    5?   4?  yes   X  X 
    8     Cash(5,2,4)       SDIRK   5    5?   4?  no    X  X 
    9     Cash(5,3,4)       SDIRK   5    5?   4?  no    X  X 
    10    Cash-Karp           ERK   6    5?   4?  no    -  - 
    11    Fehlberg            ERK   6    5?   5?  no    -  - 
    12    ARK4(3)6L[2]SA      ERK   6    4?   3   yes   -  - 
    13    ARK4(3)6L[2]SA   ESDIRK   6    4    3   yes   X  X 
    14    Dormand-Prince      ERK   7    5?   4?  no    -  - 
    15    Ismail(7,4,5)    ESDIRK   7    7?   5?  no    X  - 
    16    ARK5(4)8L[2]SA      ERK   8    5?   4?  yes   -  - 
    17    ARK5(4)8L[2]SA   ESDIRK   8    5?   4   yes   X  - 
  -----------------------------------------------------------

---------------------------------------------------------------*/
int ARKodeLoadButcherTable(int imethod, int *s, int *q, int *p, 
			   realtype (*A)[ARK_S_MAX], realtype *b, 
			   realtype *c, realtype *b2, 
			   realtype (*bd)[ARK_S_MAX]) 
{

  int i, j;

  /* initialize output tables to zero */
  for (i=0; i<ARK_S_MAX; i++) {
    b[i]  = ZERO;
    c[i]  = ZERO;
    b2[i] = ZERO;
    for (j=0; j<ARK_S_MAX; j++) {
      A[i][j]  = ZERO;
      bd[i][j] = ZERO;
    }
  }

  /* fill in coefficients based on method name */
  switch(imethod) {

  case(0):    /* Heun-Euler-ERK */
    *s = 2;
    *q = 2;
    *p = 1;
      
    A[1][0] = RCONST(1.0);

    b[0] = RCONST(0.5);
    b[1] = RCONST(0.5);

    b2[0] = RCONST(1.0);

    c[1] = RCONST(1.0);
    break;

  case(1):    /* Billington-SDIRK */
    *s = 3;
    *q = 3;
    *p = 2;

    A[0][0] = RCONST(0.292893218813);
    A[1][0] = RCONST(0.798989873223);
    A[1][1] = RCONST(0.292893218813);
    A[2][0] = RCONST(0.740789228841);
    A[2][1] = RCONST(0.259210771159);
    A[2][2] = RCONST(0.292893218813);

    b[0] = RCONST(0.691665115992);
    b[1] = RCONST(0.503597029883);
    b[2] = RCONST(-0.195262145876);

    b2[0] = RCONST(0.740789228840);
    b2[1] = RCONST(0.259210771159);

    c[0] = RCONST(0.292893218813);
    c[1] = RCONST(1.091883092037);
    c[2] = RCONST(1.292893218813);
    break;

  case(2):    /* TRBDF2-ESDIRK */
    *s = 3;
    *q = 3;
    *p = 2;

    A[1][0] = RCONST((2.0-SQRT2))/RCONST(2.0);
    A[1][1] = RCONST((2.0-SQRT2))/RCONST(2.0);
    A[2][0] = RCONST(SQRT2)/RCONST(4.0);
    A[2][1] = RCONST(SQRT2)/RCONST(4.0);
    A[2][2] = RCONST((2.0-SQRT2))/RCONST(2.0);

    b[0] = RCONST((1.0-SQRT2)/RCONST(4.0))/RCONST(3.0);
    b[1] = RCONST((3.0*SQRT2)/RCONST(4.0+1.0))/RCONST(3.0);
    b[2] = RCONST((2.0-SQRT2))/RCONST(6.0);

    b2[0] = RCONST(SQRT2)/RCONST(4.0);
    b2[1] = RCONST(SQRT2)/RCONST(4.0);
    b2[2] = RCONST((2.0-SQRT2))/RCONST(2.0);

    c[1] = RCONST(2.0-SQRT2);
    c[2] = RCONST(1.0);
    break;

  case(3):    /* TRX2-ESDIRK */
    *s = 3;
    *q = 3;
    *p = 2;

    A[1][0] = RCONST(0.25);
    A[1][1] = RCONST(0.25);
    A[2][0] = RCONST(0.25);
    A[2][1] = RCONST(0.5);
    A[2][2] = RCONST(0.25);

    b[0] = RCONST(1.0)/RCONST(6.0);
    b[1] = RCONST(2.0)/RCONST(3.0);
    b[2] = RCONST(1.0)/RCONST(6.0);

    b2[0] = RCONST(0.25);
    b2[1] = RCONST(0.5);
    b2[2] = RCONST(0.25);

    c[1] = RCONST(0.5);
    c[2] = RCONST(1.0);
    break;

  case(4):    /* Bogacki-Shampine-ERK */
    *s = 4;
    *q = 4;
    *p = 3;
			   
    A[1][0] = RCONST(1.0)/RCONST(2.0);
    A[2][1] = RCONST(3.0)/RCONST(4.0);
    A[3][0] = RCONST(2.0)/RCONST(9.0);
    A[3][1] = RCONST(1.0)/RCONST(3.0);
    A[3][2] = RCONST(4.0)/RCONST(9.0);

    b[0] = RCONST(2.0)/RCONST(9.0);
    b[1] = RCONST(1.0)/RCONST(3.0);
    b[2] = RCONST(4.0)/RCONST(9.0);

    b2[0] = RCONST(7.0)/RCONST(24.0);
    b2[1] = RCONST(1.0)/RCONST(4.0);
    b2[2] = RCONST(1.0)/RCONST(3.0);
    b2[3] = RCONST(1.0)/RCONST(8.0);

    c[1] = RCONST(1.0)/RCONST(2.0);
    c[2] = RCONST(3.0)/RCONST(4.0);
    c[3] = RCONST(1.0);
    break;

  case(5):    /* ARK3(2)4L[2]SA-ERK */
    *s = 4;
    *q = 3;
    *p = 2;

    A[1][0] = RCONST(1767732205903.0)/RCONST(2027836641118.0);
    A[2][0] = RCONST(5535828885825.0)/RCONST(10492691773637.0);
    A[2][1] = RCONST(788022342437.0)/RCONST(10882634858940.0);
    A[3][0] = RCONST(6485989280629.0)/RCONST(16251701735622.0);
    A[3][1] = RCONST(-4246266847089.0)/RCONST(9704473918619.0);
    A[3][2] = RCONST(10755448449292.0)/RCONST(10357097424841.0);

    b[0] = RCONST(1471266399579.0)/RCONST(7840856788654.0);
    b[1] = RCONST(-4482444167858.0)/RCONST(7529755066697.0);
    b[2] = RCONST(11266239266428.0)/RCONST(11593286722821.0);
    b[3] = RCONST(1767732205903.0)/RCONST(4055673282236.0);

    b2[0] = RCONST(2756255671327.0)/RCONST(12835298489170.0);
    b2[1] = RCONST(-10771552573575.0)/RCONST(22201958757719.0);
    b2[2] = RCONST(9247589265047.0)/RCONST(10645013368117.0);
    b2[3] = RCONST(2193209047091.0)/RCONST(5459859503100.0);

    c[1] = RCONST(1767732205903.0)/RCONST(2027836641118.0);
    c[2] = RCONST(3.0)/RCONST(5.0);
    c[3] = RCONST(1.0);

    bd[0][0] = RCONST(4655552711362.0)/RCONST(22874653954995.0);
    bd[0][1] = RCONST(-215264564351.0)/RCONST(13552729205753.0);
    bd[1][0] = RCONST(-18682724506714.0)/RCONST(9892148508045.0);
    bd[1][1] = RCONST(17870216137069.0)/RCONST(13817060693119.0);
    bd[2][0] = RCONST(34259539580243.0)/RCONST(13192909600954.0);
    bd[2][1] = RCONST(-28141676662227.0)/RCONST(17317692491321.0);
    bd[3][0] = RCONST(584795268549.0)/RCONST(6622622206610.0);
    bd[3][1] = RCONST(2508943948391.0)/RCONST(7218656332882.0);
    break;

  case(6):    /* ARK3(2)4L[2]SA-ESDIRK */
    *s = 4;
    *q = 3;
    *p = 2;

    A[1][0] = RCONST(1767732205903.0)/RCONST(2027836641118.0);
    A[2][0] = RCONST(5535828885825.0)/RCONST(10492691773637.0);
    A[2][1] = RCONST(788022342437.0)/RCONST(10882634858940.0);
    A[3][0] = RCONST(6485989280629.0)/RCONST(16251701735622.0);
    A[3][1] = RCONST(-4246266847089.0)/RCONST(9704473918619.0);
    A[3][2] = RCONST(10755448449292.0)/RCONST(10357097424841.0);

    b[0] = RCONST(1471266399579.0)/RCONST(7840856788654.0);
    b[1] = RCONST(-4482444167858.0)/RCONST(7529755066697.0);
    b[2] = RCONST(11266239266428.0)/RCONST(11593286722821.0);
    b[3] = RCONST(1767732205903.0)/RCONST(4055673282236.0);

    b2[0] = RCONST(2756255671327.0)/RCONST(12835298489170.0);
    b2[1] = RCONST(-10771552573575.0)/RCONST(22201958757719.0);
    b2[2] = RCONST(9247589265047.0)/RCONST(10645013368117.0);
    b2[3] = RCONST(2193209047091.0)/RCONST(5459859503100.0);

    c[1] = RCONST(1767732205903.0)/RCONST(2027836641118.0);
    c[2] = RCONST(3.0)/RCONST(5.0);
    c[3] = RCONST(1.0);

    bd[0][0] = RCONST(4655552711362.0)/RCONST(22874653954995.0);
    bd[0][1] = RCONST(-215264564351.0)/RCONST(13552729205753.0);
    bd[1][0] = RCONST(-18682724506714.0)/RCONST(9892148508045.0);
    bd[1][1] = RCONST(17870216137069.0)/RCONST(13817060693119.0);
    bd[2][0] = RCONST(34259539580243.0)/RCONST(13192909600954.0);
    bd[2][1] = RCONST(-28141676662227.0)/RCONST(17317692491321.0);
    bd[3][0] = RCONST(584795268549.0)/RCONST(6622622206610.0);
    bd[3][1] = RCONST(2508943948391.0)/RCONST(7218656332882.0);
    break;

  case(7):    /* SDIRK-5-4 */
    *s = 5;
    *q = 5;
    *p = 4;
    A[0][0] = RCONST(0.25);
    A[1][0] = RCONST(0.5);
    A[1][1] = RCONST(0.25);
    A[2][0] = RCONST(17.0)/RCONST(50.0);
    A[2][1] = RCONST(-1.0)/RCONST(25.0);
    A[2][2] = RCONST(0.25);
    A[3][0] = RCONST(371.0)/RCONST(1360.0);
    A[3][1] = RCONST(-137.0)/RCONST(2720.0);
    A[3][2] = RCONST(15.0)/RCONST(544.0);
    A[3][3] = RCONST(0.25);
    A[4][0] = RCONST(25.0)/RCONST(24.0);
    A[4][1] = RCONST(-49.0)/RCONST(48.0);
    A[4][2] = RCONST(125.0)/RCONST(16.0);
    A[4][3] = RCONST(-85.0)/RCONST(12.0);
    A[4][4] = RCONST(0.25);

    b[0] = RCONST(25.0)/RCONST(24.0);
    b[1] = RCONST(-49.0)/RCONST(48.0);
    b[2] = RCONST(125.0)/RCONST(16.0);
    b[3] = RCONST(-85.0)/RCONST(12.0);
    b[4] = RCONST(0.25);

    b2[0] = RCONST(59.0)/RCONST(48.0);
    b2[1] = RCONST(-17.0)/RCONST(96.0);
    b2[2] = RCONST(225.0)/RCONST(32.0);
    b2[3] = RCONST(-85.0)/RCONST(12.0);

    c[0] = RCONST(0.25);
    c[1] = RCONST(0.75);
    c[2] = RCONST(11.0)/RCONST(20.0);
    c[3] = RCONST(0.5);
    c[4] = RCONST(1.0);

    bd[0][0] = RCONST(11.0)/RCONST(3.0);
    bd[0][1] = RCONST(-463.0)/RCONST(72.0);
    bd[0][2] = RCONST(217.0)/RCONST(36.0);
    bd[0][3] = RCONST(-20.0)/RCONST(9.0);
    bd[1][0] = RCONST(11.0)/RCONST(2.0);
    bd[1][1] = RCONST(-385.0)/RCONST(16.0);
    bd[1][2] = RCONST(661.0)/RCONST(24.0);
    bd[1][3] = RCONST(-10.0);
    bd[2][0] = RCONST(-125.0)/RCONST(18.0);
    bd[2][1] = RCONST(20125.0)/RCONST(432.0);
    bd[2][2] = RCONST(-8875.0)/RCONST(216.0);
    bd[2][3] = RCONST(250.0)/RCONST(27.0);
    bd[3][1] = RCONST(-85.0)/RCONST(4.0);
    bd[3][2] = RCONST(85.0)/RCONST(6.0);
    bd[4][0] = RCONST(-11.0)/RCONST(9.0);
    bd[4][1] = RCONST(557.0)/RCONST(108.0);
    bd[4][2] = RCONST(-359.0)/RCONST(54.0);
    bd[4][3] = RCONST(80.0)/RCONST(27.0);
    break;

  case(8):    /* Cash(5,2,4)-SDIRK */
    *s = 5;
    *q = 5;
    *p = 4;
    A[0][0] = RCONST(0.435866521508);
    A[1][0] = RCONST(-1.13586652150);
    A[1][1] = RCONST(0.435866521508);
    A[2][0] = RCONST(1.08543330679);
    A[2][1] = RCONST(-0.721299828287);
    A[2][2] = RCONST(0.435866521508);
    A[3][0] = RCONST(0.416349501547);
    A[3][1] = RCONST(0.190984004184);
    A[3][2] = RCONST(-0.118643265417);
    A[3][3] = RCONST(0.435866521508);
    A[4][0] = RCONST(0.896869652944);
    A[4][1] = RCONST(0.0182725272734);
    A[4][2] = RCONST(-0.0845900310706);
    A[4][3] = RCONST(-0.266418670647);
    A[4][4] = RCONST(0.435866521508);

    b[0] = RCONST(0.896869652944);
    b[1] = RCONST(0.0182725272734);
    b[2] = RCONST(-0.0845900310706);
    b[3] = RCONST(-0.266418670647);
    b[4] = RCONST(0.435866521508);

    b2[0] = RCONST((-0.7-0.5))/RCONST((-0.7-0.435866521508));
    b2[1] = RCONST((0.5-0.435866521508))/RCONST((-0.7-0.435866521508));

    c[0] = RCONST(0.435866521508);
    c[1] = RCONST(-0.7);
    c[2] = RCONST(0.8);
    c[3] = RCONST(0.924556761814);
    c[4] = RCONST(1.0);
    break;

  case(9):    /* Cash(5,3,4)-SDIRK */
    *s = 5;
    *q = 5;
    *p = 4;
    A[0][0] = RCONST(0.435866521508);
    A[1][0] = RCONST(-1.13586652150);
    A[1][1] = RCONST(0.435866521508);
    A[2][0] = RCONST(1.08543330679);
    A[2][1] = RCONST(-0.721299828287);
    A[2][2] = RCONST(0.435866521508);
    A[3][0] = RCONST(0.416349501547);
    A[3][1] = RCONST(0.190984004184);
    A[3][2] = RCONST(-0.118643265417);
    A[3][3] = RCONST(0.435866521508);
    A[4][0] = RCONST(0.896869652944);
    A[4][1] = RCONST(0.0182725272734);
    A[4][2] = RCONST(-0.0845900310706);
    A[4][3] = RCONST(-0.266418670647);
    A[4][4] = RCONST(0.435866521508);

    b[0] = RCONST(0.896869652944);
    b[1] = RCONST(0.0182725272734);
    b[2] = RCONST(-0.0845900310706);
    b[3] = RCONST(-0.266418670647);
    b[4] = RCONST(0.435866521508);

    b2[0] = RCONST(0.776691932910);
    b2[1] = RCONST(0.0297472791484);
    b2[2] = RCONST(-0.0267440239074);
    b2[3] = RCONST(0.220304811849);

    c[0] = RCONST(0.435866521508);
    c[1] = RCONST(-0.7);
    c[2] = RCONST(0.8);
    c[3] = RCONST(0.924556761814);
    c[4] = RCONST(1.0);
    break;

  case(10):    /* Cash-Karp-ERK */
    *s = 6;
    *q = 5;
    *p = 4;

    A[1][0] = RCONST(1.0)/RCONST(5.0);
    A[2][0] = RCONST(3.0)/RCONST(40.0);
    A[2][1] = RCONST(9.0)/RCONST(40.0);
    A[3][0] = RCONST(3.0)/RCONST(10.0);
    A[3][1] = RCONST(-9.0)/RCONST(10.0);
    A[3][2] = RCONST(6.0)/RCONST(5.0);
    A[4][0] = RCONST(-11.0)/RCONST(54.0);
    A[4][1] = RCONST(5.0)/RCONST(2.0);
    A[4][2] = RCONST(-70.0)/RCONST(27.0);
    A[4][3] = RCONST(35.0)/RCONST(27.0);
    A[5][0] = RCONST(1631.0)/RCONST(55296.0);
    A[5][1] = RCONST(175.0)/RCONST(512.0);
    A[5][2] = RCONST(575.0)/RCONST(13824.0);
    A[5][3] = RCONST(44275.0)/RCONST(110592.0);
    A[5][4] = RCONST(253.0)/RCONST(4096.0);

    b[0] = RCONST(37.0)/RCONST(378.0);
    b[2] = RCONST(250.0)/RCONST(621.0);
    b[3] = RCONST(125.0)/RCONST(594.0);
    b[5] = RCONST(512.0)/RCONST(1771.0);

    b2[0] = RCONST(2825.0)/RCONST(27648.0);
    b2[2] = RCONST(18575.0)/RCONST(48384.0);
    b2[3] = RCONST(13525.0)/RCONST(55296.0);
    b2[4] = RCONST(277.0)/RCONST(14336.0);
    b2[5] = RCONST(1.0)/RCONST(4.0);

    c[1] = RCONST(1.0)/RCONST(5.0);
    c[2] = RCONST(3.0)/RCONST(10.0);
    c[3] = RCONST(3.0)/RCONST(5.0);
    c[4] = RCONST(1.0);
    c[5] = RCONST(7.0)/RCONST(8.0);
    break;

  case(11):    /* Fehlberg-ERK */
    *s = 6;
    *q = 5;
    *p = 5;

    A[1][0] = RCONST(1.0)/RCONST(4.0);
    A[2][0] = RCONST(3.0)/RCONST(32.0);
    A[2][1] = RCONST(9.0)/RCONST(32.0);
    A[3][0] = RCONST(1932.0)/RCONST(2197.0);
    A[3][1] = RCONST(-7200.0)/RCONST(2197.0);
    A[3][2] = RCONST(7296.0)/RCONST(2197.0);
    A[4][0] = RCONST(439.0)/RCONST(216.0);
    A[4][1] = RCONST(-8.0);
    A[4][2] = RCONST(3680.0)/RCONST(513.0);
    A[4][3] = RCONST(-845.0)/RCONST(4104.0);
    A[5][0] = RCONST(-8.0)/RCONST(27.0);
    A[5][1] = RCONST(2.0);
    A[5][2] = RCONST(-3544.0)/RCONST(2565.0);
    A[5][3] = RCONST(1859.0)/RCONST(4104.0);
    A[5][4] = RCONST(-11.0)/RCONST(40.0);

    b[0] = RCONST(25.0)/RCONST(216.0);
    b[2] = RCONST(1408.0)/RCONST(2565.0);
    b[3] = RCONST(2197.0)/RCONST(4104.0);
    b[4] = RCONST(-1.0)/RCONST(5.0);
      
    b2[0] = RCONST(16.0)/RCONST(135.0);
    b2[2] = RCONST(6656.0)/RCONST(12825.0);
    b2[3] = RCONST(28561.0)/RCONST(56430.0);
    b2[4] = RCONST(-9.0)/RCONST(50.0);
    b2[5] = RCONST(2.0)/RCONST(55.0);

    c[1] = RCONST(1.0)/RCONST(4.0);
    c[2] = RCONST(3.0)/RCONST(8.0);
    c[3] = RCONST(12.0)/RCONST(13.0);
    c[4] = RCONST(1.0);
    c[5] = RCONST(1.0)/RCONST(2.0);
    break;

  case(12):    /* ARK4(3)6L[2]SA-ERK */
    *s = 6;
    *q = 4;
    *p = 3;

    A[1][0] = RCONST(0.5);
    A[2][0] = RCONST(13861.0)/RCONST(62500.0);
    A[2][1] = RCONST(6889.0)/RCONST(62500.0);
    A[3][0] = RCONST(-116923316275.0)/RCONST(2393684061468.0);
    A[3][1] = RCONST(-2731218467317.0)/RCONST(15368042101831.0);
    A[3][2] = RCONST(9408046702089.0)/RCONST(11113171139209.0);
    A[4][0] = RCONST(-451086348788.0)/RCONST(2902428689909.0);
    A[4][1] = RCONST(-2682348792572.0)/RCONST(7519795681897.0);
    A[4][2] = RCONST(12662868775082.0)/RCONST(11960479115383.0);
    A[4][3] = RCONST(3355817975965.0)/RCONST(11060851509271.0);
    A[5][0] = RCONST(647845179188.0)/RCONST(3216320057751.0);
    A[5][1] = RCONST(73281519250.0)/RCONST(8382639484533.0);
    A[5][2] = RCONST(552539513391.0)/RCONST(3454668386233.0);
    A[5][3] = RCONST(3354512671639.0)/RCONST(8306763924573.0);
    A[5][4] = RCONST(4040.0)/RCONST(17871.0);

    b[0] = RCONST(82889.0)/RCONST(524892.0);
    b[2] = RCONST(15625.0)/RCONST(83664.0);
    b[3] = RCONST(69875.0)/RCONST(102672.0);
    b[4] = RCONST(-2260.0)/RCONST(8211.0);
    b[5] = RCONST(1.0)/RCONST(4.0);

    b2[0] = RCONST(4586570599.0)/RCONST(29645900160.0);
    b2[2] = RCONST(178811875.0)/RCONST(945068544.0);
    b2[3] = RCONST(814220225.0)/RCONST(1159782912.0);
    b2[4] = RCONST(-3700637.0)/RCONST(11593932.0);
    b2[5] = RCONST(61727.0)/RCONST(225920.0);

    c[1] = RCONST(1.0)/RCONST(2.0);
    c[2] = RCONST(83.0)/RCONST(250.0);
    c[3] = RCONST(31.0)/RCONST(50.0);
    c[4] = RCONST(17.0)/RCONST(20.0);
    c[5] = RCONST(1.0);

    bd[0][0] = RCONST(6943876665148.0)/RCONST(7220017795957.0);
    bd[0][1] = RCONST(-54480133.0)/RCONST(30881146.0);
    bd[0][2] = RCONST(6818779379841.0)/RCONST(7100303317025.0);
    bd[2][0] = RCONST(7640104374378.0)/RCONST(9702883013639.0);
    bd[2][1] = RCONST(-11436875.0)/RCONST(14766696.0);
    bd[2][2] = RCONST(2173542590792.0)/RCONST(12501825683035.0);
    bd[3][0] = RCONST(-20649996744609.0)/RCONST(7521556579894.0);
    bd[3][1] = RCONST(174696575.0)/RCONST(18121608.0);
    bd[3][2] = RCONST(-31592104683404.0)/RCONST(5083833661969.0);
    bd[4][0] = RCONST(8854892464581.0)/RCONST(2390941311638.0);
    bd[4][1] = RCONST(-12120380.0)/RCONST(966161.0);
    bd[4][2] = RCONST(61146701046299.0)/RCONST(7138195549469.0);
    bd[5][0] = RCONST(-11397109935349.0)/RCONST(6675773540249.0);
    bd[5][1] = RCONST(3843.0)/RCONST(706.0);
    bd[5][2] = RCONST(-17219254887155.0)/RCONST(4939391667607.0);
    break;

  case(13):    /* ARK4(3)6L[2]SA-ESDIRK */
    *s = 6;
    *q = 4;
    *p = 3;

    A[1][0] = RCONST(1.0)/RCONST(4.0);
    A[1][1] = RCONST(1.0)/RCONST(4.0);
    A[2][0] = RCONST(8611.0)/RCONST(62500.0);
    A[2][1] = RCONST(-1743.0)/RCONST(31250.0);
    A[2][2] = RCONST(1.0)/RCONST(4.0);
    A[3][0] = RCONST(5012029.0)/RCONST(34652500.0);
    A[3][1] = RCONST(-654441.0)/RCONST(2922500.0);
    A[3][2] = RCONST(174375.0)/RCONST(388108.0);
    A[3][3] = RCONST(1.0)/RCONST(4.0);
    A[4][0] = RCONST(15267082809.0)/RCONST(155376265600.0);
    A[4][1] = RCONST(-71443401.0)/RCONST(120774400.0);
    A[4][2] = RCONST(730878875.0)/RCONST(902184768.0);
    A[4][3] = RCONST(2285395.0)/RCONST(8070912.0);
    A[4][4] = RCONST(1.0)/RCONST(4.0);
    A[5][0] = RCONST(82889.0)/RCONST(524892.0);
    A[5][2] = RCONST(15625.0)/RCONST(83664.0);
    A[5][3] = RCONST(69875.0)/RCONST(102672.0);
    A[5][4] = RCONST(-2260.0)/RCONST(8211.0);
    A[5][5] = RCONST(1.0)/RCONST(4.0);

    b[0] = RCONST(82889.0)/RCONST(524892.0);
    b[2] = RCONST(15625.0)/RCONST(83664.0);
    b[3] = RCONST(69875.0)/RCONST(102672.0);
    b[4] = RCONST(-2260.0)/RCONST(8211.0);
    b[5] = RCONST(1.0)/RCONST(4.0);

    c[1] = RCONST(1.0)/RCONST(2.0);
    c[2] = RCONST(83.0)/RCONST(250.0);
    c[3] = RCONST(31.0)/RCONST(50.0);
    c[4] = RCONST(17.0)/RCONST(20.0);
    c[5] = RCONST(1.0);

    b2[0] = RCONST(4586570599.0)/RCONST(29645900160.0);
    b2[2] = RCONST(178811875.0)/RCONST(945068544.0);
    b2[3] = RCONST(814220225.0)/RCONST(1159782912.0);
    b2[4] = RCONST(-3700637.0)/RCONST(11593932.0);
    b2[5] = RCONST(61727.0)/RCONST(225920.0);

    bd[0][0] = RCONST(6943876665148.0)/RCONST(7220017795957.0);
    bd[0][1] = RCONST(-54480133.0)/RCONST(30881146.0);
    bd[0][2] = RCONST(6818779379841.0)/RCONST(7100303317025.0);
    bd[2][0] = RCONST(7640104374378.0)/RCONST(9702883013639.0);
    bd[2][1] = RCONST(-11436875.0)/RCONST(14766696.0);
    bd[2][2] = RCONST(2173542590792.0)/RCONST(12501825683035.0);
    bd[3][0] = RCONST(-20649996744609.0)/RCONST(7521556579894.0);
    bd[3][1] = RCONST(174696575.0)/RCONST(18121608.0);
    bd[3][2] = RCONST(-31592104683404.0)/RCONST(5083833661969.0);
    bd[4][0] = RCONST(8854892464581.0)/RCONST(2390941311638.0);
    bd[4][1] = RCONST(-12120380.0)/RCONST(966161.0);
    bd[4][2] = RCONST(61146701046299.0)/RCONST(7138195549469.0);
    bd[5][0] = RCONST(-11397109935349.0)/RCONST(6675773540249.0);
    bd[5][1] = RCONST(3843.0)/RCONST(706.0);
    bd[5][2] = RCONST(-17219254887155.0)/RCONST(4939391667607.0);
    break;

  case(14):    /* Dormand-Prince-ERK */
    *s = 7;
    *q = 5;
    *p = 4;

    A[1][0] = RCONST(1.0)/RCONST(5.0);
    A[2][0] = RCONST(3.0)/RCONST(40.0);
    A[2][1] = RCONST(9.0)/RCONST(40.0);
    A[3][0] = RCONST(44.0)/RCONST(45.0);
    A[3][1] = RCONST(-56.0)/RCONST(15.0);
    A[3][2] = RCONST(32.0)/RCONST(9.0);
    A[4][0] = RCONST(19372.0)/RCONST(6561.0);
    A[4][1] = RCONST(-25360.0)/RCONST(2187.0);
    A[4][2] = RCONST(64448.0)/RCONST(6561.0);
    A[4][3] = RCONST(-212.0)/RCONST(729.0);
    A[5][0] = RCONST(9017.0)/RCONST(3168.0);
    A[5][1] = RCONST(-355.0)/RCONST(33.0);
    A[5][2] = RCONST(46732.0)/RCONST(5247.0);
    A[5][3] = RCONST(49.0)/RCONST(176.0);
    A[5][4] = RCONST(-5103.0)/RCONST(18656.0);
    A[6][0] = RCONST(35.0)/RCONST(384.0);
    A[6][2] = RCONST(500.0)/RCONST(1113.0);
    A[6][3] = RCONST(125.0)/RCONST(192.0);
    A[6][4] = RCONST(-2187.0)/RCONST(6784.0);
    A[6][5] = RCONST(11.0)/RCONST(84.0);

    b[0] = RCONST(35.0)/RCONST(384.0);
    b[2] = RCONST(500.0)/RCONST(1113.0);
    b[3] = RCONST(125.0)/RCONST(192.0);
    b[4] = RCONST(-2187.0)/RCONST(6784.0);
    b[5] = RCONST(11.0)/RCONST(84.0);

    b2[0] = RCONST(5179.0)/RCONST(57600.0);
    b2[2] = RCONST(7571.0)/RCONST(16695.0);
    b2[3] = RCONST(393.0)/RCONST(640.0);
    b2[4] = RCONST(-92097.0)/RCONST(339200.0);
    b2[5] = RCONST(187.0)/RCONST(2100.0);
    b2[6] = RCONST(1.0)/RCONST(40.0);

    c[1] = RCONST(1.0)/RCONST(5.0);
    c[2] = RCONST(3.0)/RCONST(10.0);
    c[3] = RCONST(4.0)/RCONST(5.0);
    c[4] = RCONST(8.0)/RCONST(9.0);
    c[5] = RCONST(1.0);
    c[6] = RCONST(1.0);
    break;

  case(15):    /* Ismail(7,4,5)-ESDIRK */
    *s = 7;
    *q = 7;
    *p = 5;

    A[1][0] = RCONST(0.28589);
    A[1][1] = RCONST(0.28589);
    A[2][0] = RCONST(0.142945);
    A[2][1] = RCONST(0.924011005);
    A[2][2] = RCONST(0.28589);
    A[3][0] = RCONST(0.16803599);
    A[3][1] = RCONST(-0.049416510);
    A[3][2] = RCONST(-0.004509476);
    A[3][3] = RCONST(0.28589);
    A[4][0] = RCONST(0.182315);
    A[4][1] = RCONST(-0.112951603);
    A[4][2] = RCONST(-0.027793233);
    A[4][3] = RCONST(0.422539833);
    A[4][4] = RCONST(0.28589);
    A[5][0] = RCONST(0.24756392);
    A[5][1] = RCONST(-0.425378071);
    A[5][2] = RCONST(-0.107036282);
    A[5][3] = RCONST(0.395700134);
    A[5][4] = RCONST(0.503260302);
    A[5][5] = RCONST(0.28589);
    A[6][0] = RCONST(0.13001804);
    A[6][2] = RCONST(-0.019290177);
    A[6][3] = RCONST(0.535386266);
    A[6][4] = RCONST(0.234313169);
    A[6][5] = RCONST(-0.166317293);
    A[6][6] = RCONST(0.28589);

    b[0] = RCONST(0.13001804);
    b[2] = RCONST(-0.019290177);
    b[3] = RCONST(0.535386266);
    b[4] = RCONST(0.234313169);
    b[5] = RCONST(-0.166317293);
    b[6] = RCONST(0.28589);

    b2[0] = RCONST(0.094388663);
    b2[2] = RCONST(-0.039782614);
    b2[3] = RCONST(0.745608552);
    b2[4] = RCONST(-0.505129807);
    b2[5] = RCONST(0.704915206);

    c[1] = RCONST(0.57178);
    c[2] = RCONST(1.352846);
    c[3] = RCONST(0.4);
    c[4] = RCONST(0.75);
    c[5] = RCONST(0.9);
    c[6] = RCONST(1.0);
    break;

  case(16):    /* ARK5(4)8L[2]SA-ERK */
    *s = 8;
    *q = 5;
    *p = 4;

    A[1][0] = RCONST(41.0)/RCONST(100.0);
    A[2][0] = RCONST(367902744464.0)/RCONST(2072280473677.0);
    A[2][1] = RCONST(677623207551.0)/RCONST(8224143866563.0);
    A[3][0] = RCONST(1268023523408.0)/RCONST(10340822734521.0);
    A[3][2] = RCONST(1029933939417.0)/RCONST(13636558850479.0);
    A[4][0] = RCONST(14463281900351.0)/RCONST(6315353703477.0);
    A[4][2] = RCONST(66114435211212.0)/RCONST(5879490589093.0);
    A[4][3] = RCONST(-54053170152839.0)/RCONST(4284798021562.0);
    A[5][0] = RCONST(14090043504691.0)/RCONST(34967701212078.0);
    A[5][2] = RCONST(15191511035443.0)/RCONST(11219624916014.0);
    A[5][3] = RCONST(-18461159152457.0)/RCONST(12425892160975.0);
    A[5][4] = RCONST(-281667163811.0)/RCONST(9011619295870.0);
    A[6][0] = RCONST(19230459214898.0)/RCONST(13134317526959.0);
    A[6][2] = RCONST(21275331358303.0)/RCONST(2942455364971.0);
    A[6][3] = RCONST(-38145345988419.0)/RCONST(4862620318723.0);
    A[6][4] = RCONST(-1.0)/RCONST(8.0);
    A[6][5] = RCONST(-1.0)/RCONST(8.0);
    A[7][0] = RCONST(-19977161125411.0)/RCONST(11928030595625.0);
    A[7][2] = RCONST(-40795976796054.0)/RCONST(6384907823539.0);
    A[7][3] = RCONST(177454434618887.0)/RCONST(12078138498510.0);
    A[7][4] = RCONST(782672205425.0)/RCONST(8267701900261.0);
    A[7][5] = RCONST(-69563011059811.0)/RCONST(9646580694205.0);
    A[7][6] = RCONST(7356628210526.0)/RCONST(4942186776405.0);

    b[0] = RCONST(-872700587467.0)/RCONST(9133579230613.0);
    b[3] = RCONST(22348218063261.0)/RCONST(9555858737531.0);
    b[4] = RCONST(-1143369518992.0)/RCONST(8141816002931.0);
    b[5] = RCONST(-39379526789629.0)/RCONST(19018526304540.0);
    b[6] = RCONST(32727382324388.0)/RCONST(42900044865799.0);
    b[7] = RCONST(41.0)/RCONST(200.0);

    b2[0] = RCONST(-975461918565.0)/RCONST(9796059967033.0);
    b2[3] = RCONST(78070527104295.0)/RCONST(32432590147079.0);
    b2[4] = RCONST(-548382580838.0)/RCONST(3424219808633.0);
    b2[5] = RCONST(-33438840321285.0)/RCONST(15594753105479.0);
    b2[6] = RCONST(3629800801594.0)/RCONST(4656183773603.0);
    b2[7] = RCONST(4035322873751.0)/RCONST(18575991585200.0);

    c[1] = RCONST(41.0)/RCONST(100.0);
    c[2] = RCONST(2935347310677.0)/RCONST(11292855782101.0);
    c[3] = RCONST(1426016391358.0)/RCONST(7196633302097.0);
    c[4] = RCONST(92.0)/RCONST(100.0);
    c[5] = RCONST(24.0)/RCONST(100.0);
    c[6] = RCONST(3.0)/RCONST(5.0);
    c[7] = RCONST(1.0);

    bd[0][0] = RCONST(-17674230611817.0)/RCONST(10670229744614.0);
    bd[0][1] = RCONST(43486358583215.0)/RCONST(12773830924787.0);
    bd[0][2] = RCONST(-9257016797708.0)/RCONST(5021505065439.0);
    bd[3][0] = RCONST(65168852399939.0)/RCONST(7868540260826.0);
    bd[3][1] = RCONST(-91478233927265.0)/RCONST(11067650958493.0);
    bd[3][2] = RCONST(26096422576131.0)/RCONST(11239449250142.0);
    bd[4][0] = RCONST(15494834004392.0)/RCONST(5936557850923.0);
    bd[4][1] = RCONST(-79368583304911.0)/RCONST(10890268929626.0);
    bd[4][2] = RCONST(92396832856987.0)/RCONST(20362823103730.0);
    bd[5][0] = RCONST(-99329723586156.0)/RCONST(26959484932159.0);
    bd[5][1] = RCONST(-12239297817655.0)/RCONST(9152339842473.0);
    bd[5][2] = RCONST(30029262896817.0)/RCONST(10175596800299.0);
    bd[6][0] = RCONST(-19024464361622.0)/RCONST(5461577185407.0);
    bd[6][1] = RCONST(115839755401235.0)/RCONST(10719374521269.0);
    bd[6][2] = RCONST(-26136350496073.0)/RCONST(3983972220547.0);
    bd[7][0] = RCONST(-6511271360970.0)/RCONST(6095937251113.0);
    bd[7][1] = RCONST(5843115559534.0)/RCONST(2180450260947.0);
    bd[7][2] = RCONST(-5289405421727.0)/RCONST(3760307252460.0);
    break;

  case(17):    /* ARK5(4)8L[2]SA-ESDIRK */
    *s = 8;
    *q = 5;
    *p = 4;

    A[1][0] = RCONST(41.0)/RCONST(200.0);
    A[1][1] = RCONST(41.0)/RCONST(200.0);
    A[2][0] = RCONST(41.0)/RCONST(400.0);
    A[2][1] = RCONST(-567603406766.0)/RCONST(11931857230679.0);
    A[2][2] = RCONST(41.0)/RCONST(200.0);
    A[3][0] = RCONST(683785636431.0)/RCONST(9252920307686.0);
    A[3][2] = RCONST(-110385047103.0)/RCONST(1367015193373.0);
    A[3][3] = RCONST(41.0)/RCONST(200.0);
    A[4][0] = RCONST(3016520224154.0)/RCONST(10081342136671.0);
    A[4][2] = RCONST(30586259806659.0)/RCONST(12414158314087.0);
    A[4][3] = RCONST(-22760509404356.0)/RCONST(11113319521817.0);
    A[4][4] = RCONST(41.0)/RCONST(200.0);
    A[5][0] = RCONST(218866479029.0)/RCONST(1489978393911.0);
    A[5][2] = RCONST(638256894668.0)/RCONST(5436446318841.0);
    A[5][3] = RCONST(-1179710474555.0)/RCONST(5321154724896.0);
    A[5][4] = RCONST(-60928119172.0)/RCONST(8023461067671.0);
    A[5][5] = RCONST(41.0)/RCONST(200.0);
    A[6][0] = RCONST(1020004230633.0)/RCONST(5715676835656.0);
    A[6][2] = RCONST(25762820946817.0)/RCONST(25263940353407.0);
    A[6][3] = RCONST(-2161375909145.0)/RCONST(9755907335909.0);
    A[6][4] = RCONST(-211217309593.0)/RCONST(5846859502534.0);
    A[6][5] = RCONST(-4269925059573.0)/RCONST(7827059040749.0);
    A[6][6] = RCONST(41.0)/RCONST(200.0);
    A[7][0] = RCONST(-872700587467.0)/RCONST(9133579230613.0);
    A[7][3] = RCONST(22348218063261.0)/RCONST(9555858737531.0);
    A[7][4] = RCONST(-1143369518992.0)/RCONST(8141816002931.0);
    A[7][5] = RCONST(-39379526789629.0)/RCONST(19018526304540.0);
    A[7][6] = RCONST(32727382324388.0)/RCONST(42900044865799.0);
    A[7][7] = RCONST(41.0)/RCONST(200.0);

    b[0] = RCONST(-872700587467.0)/RCONST(9133579230613.0);
    b[3] = RCONST(22348218063261.0)/RCONST(9555858737531.0);
    b[4] = RCONST(-1143369518992.0)/RCONST(8141816002931.0);
    b[5] = RCONST(-39379526789629.0)/RCONST(19018526304540.0);
    b[6] = RCONST(32727382324388.0)/RCONST(42900044865799.0);
    b[7] = RCONST(41.0)/RCONST(200.0);

    b2[0] = RCONST(-975461918565.0)/RCONST(9796059967033.0);
    b2[3] = RCONST(78070527104295.0)/RCONST(32432590147079.0);
    b2[4] = RCONST(-548382580838.0)/RCONST(3424219808633.0);
    b2[5] = RCONST(-33438840321285.0)/RCONST(15594753105479.0);
    b2[6] = RCONST(3629800801594.0)/RCONST(4656183773603.0);
    b2[7] = RCONST(4035322873751.0)/RCONST(18575991585200.0);

    c[1] = RCONST(41.0)/RCONST(100.0);
    c[2] = RCONST(2935347310677.0)/RCONST(11292855782101.0);
    c[3] = RCONST(1426016391358.0)/RCONST(7196633302097.0);
    c[4] = RCONST(92.0)/RCONST(100.0);
    c[5] = RCONST(24.0)/RCONST(100.0);
    c[6] = RCONST(3.0)/RCONST(5.0);
    c[7] = RCONST(1.0);

    bd[0][0] = RCONST(-17674230611817.0)/RCONST(10670229744614.0);
    bd[0][1] = RCONST(43486358583215.0)/RCONST(12773830924787.0);
    bd[0][2] = RCONST(-9257016797708.0)/RCONST(5021505065439.0);
    bd[3][0] = RCONST(65168852399939.0)/RCONST(7868540260826.0);
    bd[3][1] = RCONST(-91478233927265.0)/RCONST(11067650958493.0);
    bd[3][2] = RCONST(26096422576131.0)/RCONST(11239449250142.0);
    bd[4][0] = RCONST(15494834004392.0)/RCONST(5936557850923.0);
    bd[4][1] = RCONST(-79368583304911.0)/RCONST(10890268929626.0);
    bd[4][2] = RCONST(92396832856987.0)/RCONST(20362823103730.0);
    bd[5][0] = RCONST(-99329723586156.0)/RCONST(26959484932159.0);
    bd[5][1] = RCONST(-12239297817655.0)/RCONST(9152339842473.0);
    bd[5][2] = RCONST(30029262896817.0)/RCONST(10175596800299.0);
    bd[6][0] = RCONST(-19024464361622.0)/RCONST(5461577185407.0);
    bd[6][1] = RCONST(115839755401235.0)/RCONST(10719374521269.0);
    bd[6][2] = RCONST(-26136350496073.0)/RCONST(3983972220547.0);
    bd[7][0] = RCONST(-6511271360970.0)/RCONST(6095937251113.0);
    bd[7][1] = RCONST(5843115559534.0)/RCONST(2180450260947.0);
    bd[7][2] = RCONST(-5289405421727.0)/RCONST(3760307252460.0);
    break;

  default:

    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeGetButcherTable", "Unknown Butcher table");
    return(ARK_ILL_INPUT);

  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
