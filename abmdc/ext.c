#include  <math.h>
#include  <stdlib.h>
#include  "headers/md.h"


/*********************************************************/
/* Approximated Erfc() in double precision calculation   */
/* REF: "Numerical Recipes in C" [Japanese edition]      */
/* Gijyutsu Hyouron Sya (1993), p176                     */
/* note: this function has only 10^(-7) accuracy         */
/*********************************************************/
double   erfcc(double x)
{
  double t, z, ans;

  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
        t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
        t*(-0.82215223+t*0.17087277)))))))));

  return x >= 0.0 ? ans : 2.0 - ans;

  /*
  double ans, x2;
  double Pi4;
  double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14;

  Pi4 = 12.56637061435917;
  c1 =  0.3183098861837907;
  c2 =  0.7788007830714049;
  c3 =  0.3678794411714423;
  c4 =  0.1053992245618643;
  c5 =  0.01831563888873417;
  c6 =  0.001930454136227709;
  c7 =  0.0001234098040866795;
  c8 =  4.785117392129008e-6;
  c9 =  1.125351747192591e-7;
  c10 =  1.605228055185611e-9;
  c11 =  1.388794386496401e-11;
  c12 =  7.287724095819693e-14;
  c13 =  2.319522830243569e-16;
  c14 =  4.4777324417183e-19;

  x2 = x*x;

  ans = -2.0/(-1.0 + exp(Pi4*x)) +
    (c1*x*(0.5/x2 + c2/(0.25 + x2) + c3/(1.0 + x2) +
	   c4/(2.25 + x2) + c5/(4.0 + x2) + c6/(6.25 + x2) +
	   c7/(9.0 + x2) + c8/(12.25 + x2) + c9/(16.0+ x2) +
	   c10/(20.25 + x2) + c11/(25.0 + x2) + c12/(30.25 + x2) +
	   c13/(36.0 + x2) + c14/(42.25 + x2)))/exp(x2);

  return x >= 0.0 ? ans : 2.0 - ans;  
  */

}
