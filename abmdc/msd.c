#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

/* -------------------------------------------------------------------
*   Mean Squared Displacement calculation
*
*   MSD(t) = <|r(t)-r(0)|^2> => 6 D t   [Einstein's relation]
*
*   diffusion constant D = <|r(t)-r(0)|^2>/(6t)
*
* ------------------------------------------------------------------- */
void calc_msd(void)
{
  int i;
  static short cnt=0;
  double Lx2, Ly2, Lz2;
  double dx, dy, dz, dr2;
  static double rx0[MAX_ATOMS], ry0[MAX_ATOMS], rz0[MAX_ATOMS];
  static double msdx[MAX_ATOMS], msdy[MAX_ATOMS], msdz[MAX_ATOMS];

  if(cnt==0) { /* initialization of calc_msd() at first time step */
    for(i=0; i<atom.N; i++){
      rx0[i] = cell.rx[i];
      ry0[i] = cell.ry[i];
      rz0[i] = cell.rz[i];
      msdx[i] = 0.0;
      msdy[i] = 0.0;
      msdz[i] = 0.0;
      cnt = 1;
    }
    return;  /* exit this function */
  }

  for(i=0;i<atom.ntypes;i++) {  /* initialization */
    msd.value[i] = 0.0;
    msd.number_of_ion[i] = 0;
  }
  Lx2 = cell.Lx/2.0; /* in the P=const MD, Lx, Ly, and Lz will */
  Ly2 = cell.Ly/2.0; /* change with time-steps                 */
  Lz2 = cell.Lz/2.0;

  /*------ start the actual MSD calculation in time steps >= 1 ------*/
  for(i=0; i<atom.ntypes; i++){
    msd.value[i] = 0.0;
    msd.number_of_ion[i] = 0;
  }

  for(i=0; i<atom.N; i++){
    dx = rx0[i] - cell.rx[i];  /* calculate dr = r(t-1) - r(t) */
    dy = ry0[i] - cell.ry[i];
    dz = rz0[i] - cell.rz[i];

    rx0[i] = cell.rx[i];       /* record the old positions     */
    ry0[i] = cell.ry[i];
    rz0[i] = cell.rz[i];

    /* cyclic BC effect consideration for dr = r(t-1) - r(t) */
    if( dx >  Lx2 ) dx -= cell.Lx; 
    if( dx < -Lx2 ) dx += cell.Lx;
    if( dy >  Ly2 ) dy -= cell.Ly;
    if( dy < -Ly2 ) dy += cell.Ly;
    if( dz >  Lz2 ) dz -= cell.Lz;
    if( dz < -Lz2 ) dz += cell.Lz;

    msdx[i] += dx;
    msdy[i] += dy;
    msdz[i] += dz;
    dr2 = msdx[i]*msdx[i]+msdy[i]*msdy[i]+msdz[i]*msdz[i];

    msd.value[atom.type[i]] += dr2;
    msd.number_of_ion[atom.type[i]]++;
  }

  /* evaluate an ensemble average of <MSD(t)> */
  for(i=0; i<atom.ntypes; i++){
    msd.value[i] /= (double)msd.number_of_ion[i];
  } 

}
