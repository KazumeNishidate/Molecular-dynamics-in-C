#include  <stdlib.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/* -------------------------------------------------------------------
*   Mean Squared Displacement calculation
*
*   MSD(t) = <|r(t)-r(0)|^2> => 6 D t   [Einstein's relation]
*
*   diffusion constant D = <|r(t)-r(0)|^2>/(6t)
*
*   note: this is a optional calculation. to calculate MSD, just
*   uncomment out the corresponding function-call line of "newton()"
*   in "main.c".
*
*   related groval valiables: "msd" *typedef* structure (see "md.h")
* ------------------------------------------------------------------- */

void calc_msd(void)
{
  int i, atom_kind;
  double ddx, ddy, ddz;
  double rx, ry, rz;

  if(sys.step==1){ /* initialization of calc_msd() at time step = 0 */

    for(i=0; i<sys.N; i++) {
      msd.rx_old[i] = sys.rx[i]; /* record the initial positions [t=0] */
      msd.ry_old[i] = sys.ry[i];
      msd.rz_old[i] = sys.rz[i];
    }
    return;  /* exit this function */
  }

  /*------ start the actual MSD. calculation at time step >= 1 ------*/
  msd.Lx2 = sys.Lx/2.0; /* in the P=const MD, Lx, Ly, and Lz will */
  msd.Ly2 = sys.Ly/2.0; /* change with time-steps                 */
  msd.Lz2 = sys.Lz/2.0;
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    msd.value[atom_kind] = 0.0;
    msd.number_of_ion[atom_kind] = 0;
  }

  for(i=0; i<sys.N; i++) {
    rx = sys.rx[i];
    ry = sys.ry[i];
    rz = sys.rz[i];

    ddx = msd.rx_old[i] - rx;  /* calculate dr = r(t-1) - r(t) */
    ddy = msd.ry_old[i] - ry;
    ddz = msd.rz_old[i] - rz;

    msd.rx_old[i] = rx; /* record the old positions     */
    msd.ry_old[i] = ry;
    msd.rz_old[i] = rz;

    /* cyclic BC effect consideration for dr = r(t-1) - r(t) */
    if( ddx >  msd.Lx2 ) ddx -= sys.Lx; 
    if( ddx < -msd.Lx2 ) ddx += sys.Lx;
    if( ddy >  msd.Ly2 ) ddy -= sys.Ly;
    if( ddy < -msd.Ly2 ) ddy += sys.Ly;
    if( ddz >  msd.Lz2 ) ddz -= sys.Lz;
    if( ddz < -msd.Lz2 ) ddz += sys.Lz;

    msd.dx[i] += ddx;  /* sum-up the delta_r as a function       */
    msd.dy[i] += ddy;  /* of time-step, so that we can calculate */
    msd.dz[i] += ddz;  /* the MSD over Long-time steps!          */

    /* recored the calculated MSD(t) for each kind of ions */
    msd.value[sys.ion[i]] += msd.dx[i]*msd.dx[i] +
      msd.dy[i]*msd.dy[i] + msd.dz[i]*msd.dz[i];
    msd.number_of_ion[sys.ion[i]]++; /* count ions */
  }

  /* evaluate an ensemble average of MSD. = <MSD> */
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    msd.value[atom_kind] /= (double)msd.number_of_ion[atom_kind];
  } 

}

