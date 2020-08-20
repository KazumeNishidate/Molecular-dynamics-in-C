#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

void  ewald(void)
{
  clear_fcc();    /* clear the core-core interaction force */
  ewald_I();      /* 1st term (real space)                 */
  ewald_II();     /* 2nd term (reciprocal space)           */
  ewald_III();    /* 3rd and 4th terms                     */
}

/*---------------------------------------------------*/
/* select the sphere reciprocal lattice vectors      */
/*---------------------------------------------------*/
void  init_ewald(void) /* called by main() only at once */
{
  int hx, hy, hz;
  double h2, hmax, minlat;
  int cnt=0;

  /*--- EWALD parameters -------------------------------------------------*/
  minlat = find_min3(cell.Lx, cell.Ly, cell.Lz);
  ew.alpha  = 5.6/minlat;        /* alpha parameter setting : fixed */

  printf(">>> ew.alpha = %f [/a.u]\n",ew.alpha);

  ew.hmax2  = 27.0;              /* |h^2| <= 27 : presently fixed   */
  /*----------------------------------------------------------------------*/
  ew.alpha2 = ew.alpha*ew.alpha;
  hmax  = sqrt(ew.hmax2);

  for(hx=-hmax; hx<=hmax; hx++) {
    for(hy=-hmax; hy<=hmax; hy++) {
      for(hz=-hmax; hz<=hmax; hz++) {

        if(hz==0 && hy==0 && hx==0) continue; /* without the Gamma */

        h2 = hx*hx + hy*hy + hz*hz;
        if( h2 > ew.hmax2 ) continue;
	ew.kx[cnt] = PI2*((double)hx)/cell.Lx;
	ew.ky[cnt] = PI2*((double)hy)/cell.Ly;
	ew.kz[cnt] = PI2*((double)hz)/cell.Lz;
	ew.knorm2[cnt] = ew.kx[cnt]*ew.kx[cnt] + 
	  ew.ky[cnt]*ew.ky[cnt] + ew.kz[cnt]*ew.kz[cnt];
	cnt++;
      }
    }
  }

  if(cnt > EWALD_VECTORS) {
    printf("EWALD: EWALD_VECTORS = %d\n",EWALD_VECTORS);
    printf("EWALD: generated vectors = %d\n",cnt);
    printf("EWALD: inclease the EWALD_VECTORS.\n");
    exit(0);
  }
  printf(">>> hmax = %f  ewald points = %d (without the Gamma)\n",hmax,cnt);
  ew.N = cnt;

}

void  clear_fcc(void) /* core-core force initialization */
{
  int i;

  for(i=0;i<atom.N;i++) {
    wk.fccx[i] = 0.0;
    wk.fccy[i] = 0.0;
    wk.fccz[i] = 0.0;
  }
}

void  ewald_III(void)
{
  int i, j;
  double pot1=0.0, pot2=0.0;

  for(i=0;i<atom.N;i++) { /* 3rd term */
    for(j=0;j<atom.N;j++) {
      pot1 -= atom.Zv[atom.type[i]]*atom.Zv[atom.type[j]]
	*PI/(cell.vol*ew.alpha2);
    }
  }
  for(i=0;i<atom.N;i++) { /* 4th term */
    pot2 -= (2.0*ew.alpha/SQRTPI)*atom.Zv[atom.type[i]]*atom.Zv[atom.type[i]];
  }
  E.EW += 0.5*(pot1+pot2);

}

void ewald_II(void) /* G-space */
{
  double pot=0.0, ss, cc, ti;
  double dx, dy, dz, force;
  int nk, i, j;

  for(nk=0;nk<ew.N;nk++) { /* potential */
    ss = 0.0; cc = 0.0;
    for(i=0;i<atom.N;i++) {
      ti = ew.kx[nk]*cell.rx[i] 
	 + ew.ky[nk]*cell.ry[i] 
	 + ew.kz[nk]*cell.rz[i];
      cc += atom.Zv[atom.type[i]]*cos(ti);
      ss += atom.Zv[atom.type[i]]*sin(ti);
    }
    pot += (4.0*PI/cell.vol)*(1.0/ew.knorm2[nk]) *
      exp(-ew.knorm2[nk]/(4.0*ew.alpha2)) * (cc*cc + ss*ss);
  }
  E.EW += 0.5*pot;

  for(nk=0;nk<ew.N;nk++) { /* core-core force */
    for(i=0;i<atom.N;i++) {
      for(j=0;j<atom.N;j++) {
	dx = cell.rx[j]-cell.rx[i];
	dy = cell.ry[j]-cell.ry[i];
	dz = cell.rz[j]-cell.rz[i];
	ti = ew.kx[nk]*dx+ew.ky[nk]*dy+ew.kz[nk]*dz;
	force = -(4.0*PI/cell.vol)*(1.0/ew.knorm2[nk]) *
	  exp(-ew.knorm2[nk]/(4.0*ew.alpha2)) *
	  atom.Zv[atom.type[i]]*atom.Zv[atom.type[j]]*sin(ti);
	wk.fccx[i] += force*ew.kx[nk];
	wk.fccy[i] += force*ew.ky[nk];
	wk.fccz[i] += force*ew.kz[nk];
      }
    }
  }

}

void ewald_I(void) /* real-space */
{
  int i, j, cx, cy, cz;
  double dx, dy, dz, dr, dr2;
  double force;

  E.EW = 0.0; /* initialization [potential] */

  for(cx=-2;cx<=2;cx++){ /* contribution from the replica cells */
    for(cy=-2;cy<=2;cy++){
      for(cz=-2;cz<=2;cz++){

	for(i=0;i<atom.N;i++) {	/*----- core-core force ---------*/
	  for(j=0;j<atom.N;j++) {

	    dx = cell.rx[j]-cell.rx[i]+cx*cell.Lx;
	    dy = cell.ry[j]-cell.ry[i]+cy*cell.Ly;
	    dz = cell.rz[j]-cell.rz[i]+cz*cell.Lz;
	    dr2 = dx*dx + dy*dy + dz*dz;
	    if(dr2 < 1.0e-10) continue;
	    dr = sqrt(dr2);

	    force = -atom.Zv[atom.type[i]]*atom.Zv[atom.type[j]]*
	      ( (2.0*ew.alpha/SQRTPI)*exp(-ew.alpha2*dr2)/dr 
		+ erfcc(ew.alpha*dr)/dr2 ); 
	    
	    wk.fccx[i] += force*dx/dr;
	    wk.fccy[i] += force*dy/dr;
	    wk.fccz[i] += force*dz/dr;
	    E.EW += atom.Zv[atom.type[i]]*atom.Zv[atom.type[j]]
	      *erfcc(ew.alpha*dr)/dr;
	  }
	} /*-----------------------------------------------------*/
	
      }
    }
  }
  E.EW *= 0.5;

}

double find_min3(double x, double y, double z)
{
  double min=0.0;

  if(x <= y && x <= z) min = x;
  if(y <= x && y <= z) min = y;
  if(z <= x && z <= y) min = z;
  if(min == 0.0) {printf("ERROR in find_min3\n");exit(0);}
  return min;
}
