#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/* electronic potential evaluation for occupied AOs */
void calc_Fermi_E(int target_atom)
{
  int i, j, k, order;
  int NN;
  double x;

  int local_N;

  local_N = local.num_of_nn_atoms[target_atom];
  NN = local_N*4;

  /* htb.potential = 0.0; */

  /* check_H_chebyshev(); */

  /*-------------------------------------------------------*/
  for(i=0;i<NN;i++){    /* initialization of Chebyshev series */
    for(j=0;j<NN;j++){
      if(i==j){
	fermi.T1[NN*i+j] = 1.0;        /* make up T0 = I (unit matrix) */
	fermi.mat[i*NN+j] = ch.H[1]/2.0;                 /* first term */
      } else {
	fermi.T1[NN*i+j] = 0.0;
	fermi.mat[i*NN+j] = 0.0;                         /* first term */
      }
      fermi.htb[NN*i+j] = htb.mat1[NN*i+j]/ch.Emax;  /* htb = Htb/Emax */
      fermi.T2[NN*i+j] = fermi.htb[NN*i+j];          /* -1 < htb < 1   */
      fermi.mat[i*NN+j] += ch.H[2]*fermi.htb[NN*i+j];   /* second term */
    }
  }

  for(order=3;order<ch.Npl+1;order++){ /* Chebyshev series expansion */
    for(i=0;i<NN;i++){ /* T(n+1,x) = 2.0 * x * T(n,x) - T(n-1,x) */
      for(j=0;j<4;j++){ /* NN -> 4 LOCAL orbits */
	x = 0.0;
	for(k=0;k<NN;k++){
	  x += fermi.htb[i*NN+k]*fermi.T2[k*NN+j];
	} 
	fermi.T3[i*NN+j] = 2.0*x - fermi.T1[i*NN+j]; 
      }
    }

    for(i=0;i<NN;i++){
      for(j=0;j<4;j++){ /* NN -> 4 LOCAL orbits */
	fermi.mat[i*NN+j] += ch.H[order]*fermi.T3[i*NN+j];
	fermi.T1[i*NN+j] = fermi.T2[i*NN+j]; /* T(n-1,x) -> T(n,x)   */
	fermi.T2[i*NN+j] = fermi.T3[i*NN+j]; /* T(n,x)   -> T(n+1,x) */
      }
    }
  }

  /*-------------------------------------------------------*/
  x = 0.0; /* real number of electrons = fermi.electrons */
  for(i=0;i<4;i++) {
    x += fermi.mat[NN*i + i];
  } 
  fermi.number_of_electrons = x;      /* Ne = Tr[F[htb]] */
  /*-------------------------------------------------------*/

  /* evaluate the electron energy -----------------------*/
  x = 0.0; 
  for(i=0;i<4;i++){ /* NN -> 4 LOCAL orbits */
    for(k=0;k<NN;k++){
      x += htb.mat1[i*NN+k]*fermi.mat[k*NN+i];
    } 
  }
  /*-----------------------------------------------------*/

  htb.potential[target_atom] = 2.0*x;   /* Eel = 2.0*Tr[Htb*F[htb]] */
  
/* check */
/*
  printf(" htb.potential = %f\n",htb.potential/htb.eV2E);
  printf(" electrons = %f\n",fermi.electrons);
  printf(" Tr[Fermi(htb)] = %f\n",fermi.number_of_electrons);
*/

}

/*-----------------------------------------------------------------*/

void calc_Fermi_F(int target_atom)
{
  int i, j, k, order;
  int NN;
  double x;
  double xx, yy, zz;
  int local_N;

  /* --- taka (Wednesday, June 7, 2000) --- */
  double pxx=0.0,pyy=0.0,pzz=0.0;

  local_N = local.num_of_nn_atoms[target_atom];
  NN = local_N*4;

  /*--- electric FORCE evaluation -----------------------------*/
  /*=== CHEBYSHEV ===================================*/

  for(i=0;i<NN;i++){    /* initialization of Chebyshev series */
    for(j=0;j<NN;j++){
      if(i==j){
	fermi.T1[NN*i+j] = 1.0;        /* make up T0 = I (unit matrix) */
	fermi.mat[i*NN+j] = ch.F[1]/2.0;                 /* first term */
      } else {

	fermi.T1[NN*i+j] = 0.0;
	fermi.mat[i*NN+j] = 0.0;                         /* first term */

      }
      fermi.htb[NN*i+j] = htb.mat1[NN*i+j]/ch.Emax;  /* htb = Htb/Emax */
      fermi.T2[NN*i+j] = fermi.htb[NN*i+j];          /* -1 < htb < 1   */
      fermi.mat[i*NN+j] += ch.F[2]*fermi.htb[NN*i+j];   /* second term */
    }
  }

  for(order=3;order<ch.Npl+1;order++){ /* Chebyshev series expansion */
    for(i=0;i<NN;i++){     /* T(n+1,x) = 2.0 * x * T(n,x) - T(n-1,x) */
      for(j=0;j<4;j++){  /* NN -> 4 LOCAL orbits */
	x = 0.0;
	for(k=0;k<NN;k++){
	  x += fermi.htb[i*NN+k]*fermi.T2[k*NN+j];
	} 
	fermi.T3[i*NN+j] = 2.0*x - fermi.T1[i*NN+j]; 
      }
    }

    for(i=0;i<NN;i++){
      for(j=0;j<4;j++){ /* NN -> 4 LOCAL orbits */
	fermi.mat[i*NN+j] += ch.F[order]*fermi.T3[i*NN+j];
	fermi.T1[i*NN+j] = fermi.T2[i*NN+j]; /* T(n-1,x) -> T(n,x)   */
	fermi.T2[i*NN+j] = fermi.T3[i*NN+j]; /* T(n,x)   -> T(n+1,x) */
      }
    }
  }  /*=================================================*/

  /*---------------------------------------------------------------------*/
  /* Force matrix "htb.mat2" and Fermi matrix "fermi.mat" are the        */
  /* symmetric matrix. We use here                                       */
  /*       2.0*htb.mat2x[i*NN+k]*fermi.mat[k*NN+i];                      */
  /* rather than                                                         */
  /*       2.0*fermi.mat[i*NN+k]*htb.mat2x[k*NN+i];                      */
  /*                                                                     */
  /* [note] Non-zero part of the Force matrix "htb.mat2" is in the       */
  /* region of {i => 4 && j < 4} or {i < 4 && j => 4} in "htb.mat2[i,j]",*/
  /* satisfying the relation "htb.mat2[i,j] = htb.mat2[j,i]".            */
  /*                                                                     */
  /* [proof]                                                             */
  /* In[1]:= mat1={{a,b,c,d},{b,0,0,0},{c,0,0,0},{d,0,0,0}};             */
  /* In[2]:= mat2={{A,B,C,D},{B,E,F,G},{C,F,H,I},{D,G,I,J}};             */
  /* In[4]:= Sum[(mat1.mat2)[[i,i]],{i,4}]                               */
  /* Out[4]= a A + 2 b B + 2 c C + 2 d D                                 */
  /* In[5]:= 2*Apply[Plus,Table[mat1[[i,k]]*mat2[[k,i]],{k,4},{i,1}]]    */
  /* Out[5]= {2 (a A + b B + c C + d D)}                                 */
  /*                                                                     */
  /*    where "a A" is always zero.                                      */
  /*---------------------------------------------------------------------*/

  xx = 0.0;
  yy = 0.0;
  zz = 0.0;

  /* --- taka (Wednesday, June 7, 2000) --- */
  pxx = 0.0;
  pyy = 0.0;
  pzz = 0.0;

  for(i=0;i<4;i++){ /* NN -> 4 LOCAL orbits */
    for(k=4;k<NN;k++){ /* {i<4 && k<4} elements of htb.mat2 are all zero */
      xx += 2.0*htb.mat2x[i*NN+k]*fermi.mat[k*NN+i];
      yy += 2.0*htb.mat2y[i*NN+k]*fermi.mat[k*NN+i];
      zz += 2.0*htb.mat2z[i*NN+k]*fermi.mat[k*NN+i];

      /* --- taka (Wednesday, June 7, 2000) --- */
      pxx += 2.0*press.matX[i*NN+k]*fermi.mat[k*NN+i];
      pyy += 2.0*press.matY[i*NN+k]*fermi.mat[k*NN+i];
      pzz += 2.0*press.matZ[i*NN+k]*fermi.mat[k*NN+i];
    } 
  }

  xx *= -2.0;    yy *= -2.0;    zz *= -2.0;

  /*
  printf( "pxx = %f pyy = %f pzz = %f\n", pxx, pyy, pzz );
  getchar();
  */

  /* --- taka (Wednesday, June 7, 2000) --- */
  pxx *= -2.0;   pyy *= -2.0;   pzz *= -2.0;

/*
  printf("electrical F\n");
  printf("atom = %d  %f %f %f\n",target_atom, xx,yy,zz);
  printf("repulsive F\n");
  printf("atom = %d  %f %f %f\n",target_atom, sys.fx[target_atom],
	 sys.fy[target_atom],sys.fz[target_atom]);
*/

  sys.fx[target_atom] += xx;
  sys.fy[target_atom] += yy;
  sys.fz[target_atom] += zz;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.virX += pxx;
  press.virY += pyy;
  press.virZ += pzz;

}


void check_H_chebyshev(void){
  double x;
  double xx[100];
  double y;
  int i, j;

  x = -0.999999;

  for(j=0;j<199;j++){ 

    xx[0] = 1.0;
    xx[1] = x;

    for(i=2;i<ch.Npl;i++) {
      xx[i] = 2.0*x*xx[i-1] - xx[i-2];
    }

    y = ch.H[1]/2.0;

    for(i=2;i<ch.Npl+1;i++) {
      y += ch.H[i]*xx[i-1];
    }
    printf("%f %f\n",x,y);

    x += 0.01;

  }
  exit(0);

}

/*-----------------------------------------------------------------*/
