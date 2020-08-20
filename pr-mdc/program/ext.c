#include  <math.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"
#ifndef  PI
#define  PI  3.1415926535897932385
#endif

/***********************************************/ 
/* Normal random numbers by Box-Muller method  */
/* The distribution is defined by              */
/*                1                            */
/*   p(y)dy = ---------- Exp[-y^2/2] dy        */
/*             Sqrt[2Pi]                       */
/*                                             */
/*  mean value => temp                         */
/*  dispersion => bunsan                       */
/***********************************************/ 

double   nrand(double temp, double bunsan)
{
  static unsigned int seed;
  static short cnt = 0;
  static int flag_for_recall = 0;
  static double coefficient, theta;
  double nrnd_num;
  double sqrt_velocity2_per_m;

  /* set a SEED for pseudo random number sequence */
  if(cnt == 0){
    seed = (unsigned int)time( 0 );
    srand( seed );
    cnt = 1;

    /********************************************************
     * if you want to perform several MD calculations using
     * the SAME initial velocity distribution, you should
     * set a FIXED random number SEED xxx (xxx >= 0 ).
     *   # (example) to use the FIXED SEED (=100000)
     *     srand(100000);
     ********************************************************/
  }

  /* pseudo random number sequence [0,1] */
  if(flag_for_recall) {
    flag_for_recall = 0;
    nrnd_num = bunsan * coefficient * sin( theta );
    sqrt_velocity2_per_m = sqrt(sys.kB * (temp + nrnd_num)); 
    if( (rand()/(RAND_MAX + 1.0)) < 0.5){
      return sqrt_velocity2_per_m;
    } else {
      return -sqrt_velocity2_per_m;
    }
  }

  flag_for_recall = 1;
  coefficient = sqrt(-2.0 * log(1.0 - (rand() / (RAND_MAX + 1.0))));
  theta = 2.0 * PI * (rand() / (RAND_MAX + 1.0));
  nrnd_num = bunsan * coefficient * cos( theta );
  sqrt_velocity2_per_m = sqrt(sys.kB * (temp + nrnd_num));

  if( (rand()/(RAND_MAX + 1.0)) < 0.5){
    return sqrt_velocity2_per_m;
  } else {
    return -sqrt_velocity2_per_m;
  }
}

/*********************************************************/
/* Approximated Erfc() in double precision calculation   */
/* note: the valid x value range is about 0.1 < x < 6    */
/*********************************************************/
double   erfcc(double x)
{
  static double Pi4 = 12.56637061435917;  /* 4.0 x Pi */
  static double  c1 =  0.3183098861837907;
  static double  c2 =  0.7788007830714049;
  static double  c3 =  0.3678794411714423;
  static double  c4 =  0.1053992245618643;
  static double  c5 =  0.01831563888873417;
  static double  c6 =  0.001930454136227709;
  static double  c7 =  0.0001234098040866795;
  static double  c8 =  4.785117392129008e-6;
  static double  c9 =  1.125351747192591e-7;
  static double c10 =  1.605228055185611e-9;
  static double c11 =  1.388794386496401e-11;
  static double c12 =  7.287724095819693e-14;
  static double c13 =  2.319522830243569e-16;
  static double c14 =  4.4777324417183e-19;

  double x2, erfc_value;

  x2 = x*x;

  erfc_value = -2.0/(-1.0 + exp(Pi4*x)) + 
    (c1*x*(0.5/x2 + c2/(0.25 + x2) + c3/(1.0 + x2) + 
	    c4/(2.25 + x2) + c5/(4.0 + x2) + c6/(6.25 + x2) + 
	    c7/(9.0 + x2) + c8/(12.25 + x2) + c9/(16.0+ x2) + 
	    c10/(20.25 + x2) + c11/(25.0 + x2) + c12/(30.25 + x2) + 
	    c13/(36.0 + x2) + c14/(42.25 + x2)))/exp(x2);

  return x >= 0.0 ? erfc_value : 2.0 - erfc_value; 

}
