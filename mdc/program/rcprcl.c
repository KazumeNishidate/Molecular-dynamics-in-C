#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
*  calculation of the EWALD 2nd term in the reciprocal 
*  lattice space.
*****/
void	reciprocal_space(void) 
{
  short i, j;
  short h2, hx, hy, hz;
  double hh2, hhx, hhy, hhz;
  double hhx1, hhy1, hhz1;
  double ci, si, ti;
  double ex, foc;
  double v2, vp, vv;
  double pot = 0.0 ;
  double vp_ex_per_hh2, v2_ex_per_hh2;

  if(sys.a1 == 0.0 || sys.a1 > 10.0) return;   /* SKIP the EWALD */

  vv = sys.Lx*sys.Ly*sys.Lz;
  v2 = 2.0/vv;
  vp = 1.0/(PI*vv); 

  /*------------------------------------------------------------------*/
  /* semi-sphere reciprocal lattice space sum : since V_(q) = V_(-q)  */
  /*------------------------------------------------------------------*/
  for(hx=-sys.hm_sqrt; hx<=sys.hm_sqrt; hx++) { 
    for(hy=-sys.hm_sqrt; hy<=sys.hm_sqrt; hy++) {
      for(hz=0; hz<=sys.hm_sqrt; hz++) {  
        if(hz==0 && (hy<0 || (hy==0 && hx<=0))) continue;

	h2 = hx*hx + hy*hy + hz*hz;
	if( h2 == 0 || h2 > sys.hm ) continue;
	hhx = ((double)hx)/sys.Lx;
	hhy = ((double)hy)/sys.Ly;
	hhz = ((double)hz)/sys.Lz;
	hh2 = hhx*hhx + hhy*hhy + hhz*hhz;

	hhx1 = PI2*hhx;
	hhy1 = PI2*hhy;
	hhz1 = PI2*hhz;

	ex = exp( -hh2*sys.p2a2);

	/* 
	 *  For a simulation of fixed "MD cell size" system, you
	 *  can prepare a look-up table of the function 
	 *  "ex=exp( -hh2*sys.p2a2);" to get a calculational efficiency.
	 */

	vp_ex_per_hh2 = vp*(ex/hh2);
	v2_ex_per_hh2 = v2*(ex/hh2);

	sys.rcsi[sys.N] = 0.0;
	sys.rsni[sys.N] = 0.0;  

	for(j=0; j<sys.N; j++) {

	  ti = hhx1*sys.rx[j] + hhy1*sys.ry[j]+hhz1*sys.rz[j];

	  ci = ion_z(j) * cos(ti);
	  si = ion_z(j) * sin(ti);

	  sys.rcsi[j] = ci;
	  sys.rsni[j] = si;

	  sys.rcsi[sys.N] += ci;
	  sys.rsni[sys.N] += si;
	}

	/* reciprocal space sum = 2*(semi-sphere reciprocal space sum)  */
	sys.rcsi[sys.N] *= 2.0;
	sys.rsni[sys.N] *= 2.0;

	pot += vp_ex_per_hh2*
	  (sys.rcsi[sys.N]*sys.rcsi[sys.N] + 
	   sys.rsni[sys.N]*sys.rsni[sys.N]);

	for(i=0; i<sys.N; i++) {
	  foc = -v2_ex_per_hh2*    /* BUG FIX : "+" to "-" */
	    (sys.rcsi[i]*sys.rsni[sys.N] - 
	     sys.rsni[i]*sys.rcsi[sys.N]);
	  sys.fx[i] += foc * hhx;
	  sys.fy[i] += foc * hhy;
	  sys.fz[i] += foc * hhz;

	}
      }
    }
  }
  sys.pot += pot;
}

/*****
*  EWALD potential 3rd term calculation.
*
*  In the following, the 3rd term will be evaluated only at first MD
*  calculational step for sys.N(number of particles)=const. system. 
*  It should be noticed that we must calculate this value at every
*  update time-steps when we use a non-equilibrium system such as
*  grand-canonical ensemble.
*****/ 
void	reciprocal_space3(void)
{
  static double	Z2, ewald_pot3 = 0.0;
  static int cnt = 0;
  static int old_num_of_particles;
  short i;

  if(sys.a1 == 0.0 || sys.a1 > 10.0) return;  /* SKIP the EWALD */

  /*- for a non-equilibrium system [sys.N != const.] -*/
  if(sys.step == 1){
    old_num_of_particles = sys.N;
  }
  if(old_num_of_particles != sys.N) cnt = 0;
  /*--------------------------------------------------*/

  if(cnt==0){
    Z2 = 0.0;
    for(i=0; i<sys.N; i++) {
      Z2 += ion_z(i) * ion_z(i);
    }
    ewald_pot3 = -Z2*sys.a3;
    cnt = 1;
  }
  sys.pot += ewald_pot3; 
}
