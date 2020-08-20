#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
*  calculation of the EWALD 1st term
*  no replica cell consideration
*****/
void   ewald1(void) {
  int i, j;
  double rxi, ryi, rzi;
  double dx, dy, dz;
  double ZZ, dr, adr, fij;
  double erfc_adr_per_dr;

  for(i=0; i<sys.N; i++) {
    rxi = cell.rx[i];
    ryi = cell.ry[i];
    rzi = cell.rz[i];

    for(j=i+1; j<sys.N; j++) {
      dx = rxi - cell.rx[j];
      dy = ryi - cell.ry[j];
      dz = rzi - cell.rz[j];

      cyclic_bc_dr(&dx, &dy, &dz);

      dr = dx*dx + dy*dy + dz*dz; /* dr' = dr^2 */

      if(dr > sys.cutoff2) continue; 
      dr = sqrt(dr);

      ZZ = ion_z(i) * ion_z(j);

      adr = sys.a1 * dr;  
      erfc_adr_per_dr = erfcc(adr)/dr;

      sys.pot += ZZ * erfc_adr_per_dr;

      fij = ZZ* (2.0 * (sys.a1/SQRTPI) * exp(- adr * adr) +
		 erfc_adr_per_dr )/(dr*dr);

      cell.fx[i] += fij * dx;
      cell.fy[i] += fij * dy;
      cell.fz[i] += fij * dz;
      cell.fx[j] -= fij * dx;
      cell.fy[j] -= fij * dy;
      cell.fz[j] -= fij * dz;
      /* (note): fij' => -Fij/dr, fijx => -fij'*dx => -Fij*dx/dr  */

      /* VIRIAL : EWALD 1st term */
      pt.virial[0][0] += fij * dx * dx; 
      pt.virial[0][1] += fij * dx * dy; 
      pt.virial[0][2] += fij * dx * dz; 

      pt.virial[1][1] += fij * dy * dy; 
      pt.virial[1][2] += fij * dy * dz; 

      pt.virial[2][2] += fij * dz * dz; 

    }
  }
  pt.virial[1][0] = pt.virial[0][1];
  pt.virial[2][0] = pt.virial[0][2];
  pt.virial[2][1] = pt.virial[1][2];
}

/*****
*  calculation of the EWALD 2nd term in the reciprocal 
*  lattice space.
*  no replica cell consideration
*****/
void	ewald2(void) 
{
  int i, j;
  int h2, hx, hy, hz;
  double hh2, hhx, hhy, hhz;
  double hhx1, hhy1, hhz1;
  double ci, si, ti;
  double ex, foc;
  double v2, vp;
  double pot = 0.0 ;
  double vp_ex_per_hh2, v2_ex_per_hh2;
  double **g;

  g = cell.g; /* {g} matrix */

  v2 = 2.0/cell.vol;
  vp = 1.0/(PI*cell.vol); 

  /*------------------------------------------------------------------*/
  /* semi-sphere reciprocal lattice space sum : since V_(q) = V_(-q)  */
  /*------------------------------------------------------------------*/
  for(hx=-sys.hm_sqrt; hx<=sys.hm_sqrt; hx++) { 
    for(hy=-sys.hm_sqrt; hy<=sys.hm_sqrt; hy++) {
      for(hz=0; hz<=sys.hm_sqrt; hz++) {  
        if(hz==0 && (hy<0 || (hy==0 && hx<=0))) continue;

	h2 = hx*hx + hy*hy + hz*hz;
	if( h2 == 0 || h2 > sys.hm ) continue;

	/* {g} ={ b x c, c x a, a x b }/(cell volume) */
	hhx = hx*g[0][0] + hy*g[1][0] + hz*g[2][0];
	hhy = hx*g[0][1] + hy*g[1][1] + hz*g[2][1];
	hhz = hx*g[0][2] + hy*g[1][2] + hz*g[2][2];

	hh2 = hhx*hhx + hhy*hhy + hhz*hhz;

	hhx1 = PI2*hhx;
	hhy1 = PI2*hhy;
	hhz1 = PI2*hhz;

	ex = exp( -hh2 * (PIPI/(sys.a1*sys.a1)) );

	vp_ex_per_hh2 = vp*(ex/hh2);
	v2_ex_per_hh2 = v2*(ex/hh2);

	sys.rcsi[sys.N] = 0.0;
	sys.rsni[sys.N] = 0.0;  

	for(j=0; j<sys.N; j++) {

	  ti = hhx1*cell.rx[j] + hhy1*cell.ry[j] + hhz1*cell.rz[j];

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
	  foc = - v2_ex_per_hh2*
	    (sys.rcsi[i]*sys.rsni[sys.N] - 
	     sys.rsni[i]*sys.rcsi[sys.N]);
	  cell.fx[i] += foc * hhx;
	  cell.fy[i] += foc * hhy;
	  cell.fz[i] += foc * hhz;

	}
      }
    }
  }
  sys.pot += pot;

}

/*  EWALD potential 3rd term calculation. */
void	ewald3(void)
{
  static double	Z2, ewald_pot3 = 0.0;
  int i;

  Z2 = 0.0;
  for(i=0; i<sys.N; i++) {
    Z2 += ion_z(i) * ion_z(i);
  }
  ewald_pot3 = -Z2 * (sys.a1/SQRTPI);
  sys.pot += ewald_pot3;
}

/***
 *   estimation of alpha value for EWALD calculation  [optional]
***/
void calc_alpha(void)
{
  int i;

  sys.a1 = 0.01; /* initial value for alpha parameter */

  for(i=0;i<200;i++) {
    sys.pot = 0.0;
    ewald1();       /* EWALD 1st term */
    ewald2();       /*       2nd term */
    ewald3();       /*       3rd term */

    printf("%f  %f\n",sys.a1, sys.pot);
    sys.a1  += 0.01;
  }
  exit(0);
}
