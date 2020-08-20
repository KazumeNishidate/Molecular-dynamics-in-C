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
/*    srand( 100000 );  */
    cnt = 1;

    /********************************************************
     *    if you want perform several MD calculations using
     *    the SAME initial velocity distribution, you should
     *    set a FIXED random number SEED using "srand()".
     *      example) FIXED SEED (=100000)
     *      srand(100000);
     ********************************************************/
  }

/* pseudo random number sequence [0,1] */

  if(flag_for_recall) 
    {
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
