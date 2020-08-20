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
*   note: this is a optional calculation. 
*   related groval valiables: "msd" *typedef* structure (see "md.h")
* ------------------------------------------------------------------- */

void calc_msd(void)
{
  int i, atom_kind;
  double ddx, ddy, ddz;
  double rx, ry, rz;

  if(sys.step==1){ /* initialization of calc_msd() at time step = 0 */

    for(i=0; i<sys.N; i++) {
      msd.rx_old[i] = cell.rx[i]; /* record the initial positions [t=0] */
      msd.ry_old[i] = cell.ry[i];
      msd.rz_old[i] = cell.rz[i];
    }
    return;  /* exit this function */
  }

  /*------ start the actual MSD. calculation at time step >= 1 ------*/
  for(atom_kind=0; atom_kind<ctl.kinds_of_ions; atom_kind++){
    msd.value[atom_kind] = 0.0;
    msd.number_of_ion[atom_kind] = 0;
  }

  for(i=0; i<sys.N; i++) {
    rx = cell.rx[i];
    ry = cell.ry[i];
    rz = cell.rz[i];

    ddx = msd.rx_old[i] - rx;  /* calculate dr = r(t-1) - r(t) */
    ddy = msd.ry_old[i] - ry;
    ddz = msd.rz_old[i] - rz;

    msd.rx_old[i] = rx;        /* record the old positions     */
    msd.ry_old[i] = ry;
    msd.rz_old[i] = rz;

    cyclic_bc_dr(&ddx, &ddy, &ddz);     /* cyclic BC effect    */

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

