#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "headers/md.h"

void   check_bhs_dV_in_real_space(int type, int eru)
{                   /* "eru" -> angular orbital number                     */
  int dr, i;        /* "type" = type of atom. This is not an atomic number */
  double r, r2, delta_v;
  double aa[6], alpha[3];
  double core_v=0.0;

  /*--- preperation -------------------------------------------------------*/
  for(i=0;i<6;i++) {
    aa[i] =  bhs.A[type][eru][i+1];
  }
  for(i=0;i<3;i++) {
    alpha[i] = bhs.alpha[type][eru][i+1];
  }
  /*-----------------------------------------------------------------------*/

  /*--- calculate pseudo potential v.s. length in real space --------------*/
  for(dr=0;dr<26;dr++) {  /* length [atomic unit] 0.001 up to 2.5 [a.u] */
    r = (double)dr/10.0; 
    if(r == 0.0) r = 0.001; /* to avoid the under-flow at r = 0.0 */
    r2 = r*r;

    core_v = -(atom.Zv[type]/r)*
      ( bhs.core_orbit_c[type][0]* 
	ERFF( sqrt(bhs.core_orbit_alpha[type][0]) *r ) + 
        bhs.core_orbit_c[type][1]* 
	ERFF( sqrt(bhs.core_orbit_alpha[type][1]) *r ) ) ;

    delta_v = core_v;  /* pseudo potential core part */

    for(i=0;i<3;i++) {
      delta_v += (aa[i]+r2*aa[i+3])*exp(-alpha[i]*r2);
    }
    printf(" %f  %f\n",r,delta_v);  
  }                          /* numerical table for the "r v.s. potential" */
  /*-----------------------------------------------------------------------*/
}

/************************************************************************
 *  BHS cofficient transformation program "c2a"
 *  - Ref:Computer Physics Reports Vol.9 (1989) pp. 115-198,
 *    Warren E. Pickett
 *  - The program "c2a" was originally written for the FORTRAN compiler.
 ************************************************************************/
void   init_bhs_c2a(int type)
{
  int i, j, eru;
  double alpha_sum, dumy;
  double S[6][6], Q[6][6], alpha[6], C[6], A[6];

  for(eru=0;eru<atom.norbit[type];eru++){ /* eru => angular orbital number */

    /*----------- construct overlap matrix ----------------------*/
    for(i=0;i<3;i++){ /* alpha{1,2,3,4,5,6} = alpha{1,2,3,1,2,3} */
      alpha[i] = bhs.alpha[type][eru][i+1]; /* alpha is cyclic */
      alpha[i+3] = alpha[i];
    }
    for(i=0;i<6;i++) {
      A[i] = 0.0;  /* initialize */
      C[i] = bhs.c[type][eru][i];

      for(j=0;j<6;j++) {
	S[i][j] = 0.0;  /* initialize */
	Q[i][j] = 0.0;

	alpha_sum = alpha[i]+alpha[j];
	dumy = 0.25*SQRTPI/(alpha_sum*sqrt(alpha_sum));

	if(i<3&&j<3) S[i][j] = dumy;
	if(i>2&&j>2) S[i][j] = dumy*3.75/(alpha_sum*alpha_sum);
	if((i<3&&j>2) || (i>2&&j<3)) S[i][j] = dumy*1.5/alpha_sum;

      } /* end j-loop */
    } /* end i-loop */
    /*-----------------------------------------------------------*/

    /*----------- construct transformation matrix ---------------*/
    Q[0][0] = sqrt(S[0][0]); /* 1st step */

    for(j=1;j<6;j++) {
      Q[0][j] = S[0][j]/Q[0][0];
    }

    Q[1][1] = sqrt( S[1][1] - (Q[0][1]*Q[0][1]) );
    Q[1][2] = (S[1][2] - (Q[0][1]*Q[0][2]) )/Q[1][1];
    Q[2][2] = sqrt( S[2][2] - (Q[0][2]*Q[0][2]+Q[1][2]*Q[1][2]) );
    Q[1][3] = ( S[1][3] - (Q[0][1]*Q[0][3]) )/Q[1][1];
    Q[2][3] = ( S[2][3] - (Q[0][2]*Q[0][3]+Q[1][2]*Q[1][3]) )/Q[2][2];
    Q[3][3] = sqrt( S[3][3] - 
		   (Q[0][3]*Q[0][3]+Q[1][3]*Q[1][3]+Q[2][3]*Q[2][3]) );
    Q[1][4] = ( S[1][4] - Q[0][1]*Q[0][4])/Q[1][1];
    Q[2][4] = ( S[2][4] - (Q[0][2]*Q[0][4]+Q[1][2]*Q[1][4]) )/Q[2][2];
    Q[3][4] = ( S[3][4] - 
	       (Q[0][3]*Q[0][4]+Q[1][3]*Q[1][4]+Q[2][3]*Q[2][4]) )/Q[3][3];
    Q[4][4] = sqrt( S[4][4] - (Q[0][4]*Q[0][4]+Q[1][4]*Q[1][4]+
			       Q[2][4]*Q[2][4]+Q[3][4]*Q[3][4]) );
    Q[1][5] = ( S[1][5] - (Q[0][1]*Q[0][5]) )/Q[1][1];
    Q[2][5] = ( S[2][5] - 
	       (Q[0][2]*Q[0][5]+Q[1][2]*Q[1][5]) )/Q[2][2];
    Q[3][5] = ( S[3][5] - 
	       (Q[0][3]*Q[0][5]+Q[1][3]*Q[1][5]+Q[2][3]*Q[2][5]) )/Q[3][3];
    Q[4][5] = ( S[4][5] - 
	       (Q[0][4]*Q[0][5]+Q[1][4]*Q[1][5]+
		Q[2][4]*Q[2][5]+Q[3][4]*Q[3][5]) )/Q[4][4];
    Q[5][5] = sqrt( S[5][5] -
		   (Q[0][5]*Q[0][5]+Q[1][5]*Q[1][5]+
		    Q[2][5]*Q[2][5]+Q[3][5]*Q[3][5]+Q[4][5]*Q[4][5]) );

    /*-----------------------------------------------------------*/

    /*------ construct A matrix by back-transforming ------------*/
    A[5] = -C[5]/Q[5][5];
    A[4] = -(C[4]+Q[4][5]*A[5])/Q[4][4];
    A[3] = -(C[3]+Q[3][4]*A[4]+Q[3][5]*A[5])/Q[3][3];
    A[2] = -(C[2]+Q[2][3]*A[3]+Q[2][4]*A[4]+Q[2][5]*A[5])/Q[2][2];
    A[1] = -(C[1]+Q[1][2]*A[2]+Q[1][3]*A[3]+
	     Q[1][4]*A[4]+Q[1][5]*A[5])/Q[1][1];
    A[0] = -(C[0]+Q[0][1]*A[1]+Q[0][2]*A[2]+Q[0][3]*A[3]+Q[0][4]*A[4]+
	     Q[0][5]*A[5])/Q[0][0];

    for(i=0;i<6;i++) {
      bhs.A[type][eru][i+1] = A[i];
    }
    /*-----------------------------------------------------------*/
  } /* end eru-loop */

}
