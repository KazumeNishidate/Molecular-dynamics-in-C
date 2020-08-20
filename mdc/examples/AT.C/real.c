#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*-----------------------------------------------------------------------*/
/*  The function for potential and force calculation with Abell Tersoff  */ 
/*  potentail set.                                                       */
/*-----------------------------------------------------------------------*/
void	real_space(void) {
  int i, j, k;
  double c2, d2;
  double rdx, rdy, rdz;
  double rij, xij, yij, zij;
  double rik, xik, yik, zik;
  double fc_rij, drij_fc_rij;
  double fc_rik, drik_fc_rik;
  double xxx, ggg;
  double cs, si;
  double beta_n, xxx_n, xxx_n1;
  double bij, drik_bij, dtheta_bij;
  double f1, f2, f3;
  double tmp;
  double dxj_theta, dyj_theta, dzj_theta;
  double dxk_theta, dyk_theta, dzk_theta;

  rdx = sys.Lx*0.5;
  rdy = sys.Ly*0.5;
  rdz = sys.Lz*0.5;

  c2 = at.c * at.c;
  d2 = at.d * at.d;

  beta_n = exp(at.n*log(at.beta));    /* beta^n */

  /* initial set */
  sys.pot = 0.0;
  sys.virX = 0.0;
  sys.virY = 0.0;
  sys.virZ = 0.0;

  for(i=0;i<sys.N;i++) {
    for(j=0;j<sys.N;j++) {
      if(i==j) continue;
      xij = sys.rx[i] - sys.rx[j];
      yij = sys.ry[i] - sys.ry[j];
      zij = sys.rz[i] - sys.rz[j];

      /* cyclic boundary condition */
      if(xij > rdx) xij -= sys.Lx;
      if(yij > rdy) yij -= sys.Ly;
      if(zij > rdz) zij -= sys.Lz;
      if(xij < -rdx) xij += sys.Lx;
      if(yij < -rdy) yij += sys.Ly;
      if(zij < -rdz) zij += sys.Lz;

      if(xij > at.S || xij < -at.S) continue;
      if(yij > at.S || yij < -at.S) continue;
      if(zij > at.S || zij < -at.S) continue;

      rij = sqrt(xij*xij + yij*yij + zij*zij);

      /*-----------  f_c(rij) --------------------------------------*/
      fc_rij = 0.0;
      drij_fc_rij = 0.0;

      if(rij > at.S) continue;
      if(rij < at.R) fc_rij = 1.0;
      if(rij >= at.R && rij <= at.S)
	fc_rij = 0.5 + 0.5 * cos((PI*(rij - at.R))/(at.S - at.R));

      /*-----------  (d/drij)f_c(rij) ------------------------------*/
      if(rij < at.R || rij > at.S) drij_fc_rij=0.0;
      if(rij >= at.R && rij <= at.S)
	drij_fc_rij = -0.5*(PI/(at.S-at.R))*sin(PI*(rij-at.R)/(at.S-at.R));

      xxx = 0.0;
      for(k=0;k<sys.N;k++) {

	if(k==i||k==j) continue;

	xik = sys.rx[i] - sys.rx[k];
	yik = sys.ry[i] - sys.ry[k];
	zik = sys.rz[i] - sys.rz[k];

	/* cyclic boundary condition */
	if(xik > rdx) xik -= sys.Lx;
	if(yik > rdy) yik -= sys.Ly;
	if(zik > rdz) zik -= sys.Lz;
	if(xik < -rdx) xik += sys.Lx;
	if(yik < -rdy) yik += sys.Ly;
	if(zik < -rdz) zik += sys.Lz;

	if(xik > at.S || xik < -at.S) continue;
	if(yik > at.S || yik < -at.S) continue;
	if(zik > at.S || zik < -at.S) continue;

	rik = sqrt(xik*xik + yik*yik + zik*zik);

	/*-----------  f_c(rik) --------------------------------------*/
	fc_rik = 0.0;
	drik_fc_rik = 0.0;

	if(rik > at.S) continue;
	if(rik < at.R) fc_rik = 1.0;
	if(rik >= at.R && rik <= at.S)
	  fc_rik = 0.5 + 0.5 * cos((PI*(rik - at.R))/(at.S - at.R));

	/*-----------  (d/drik)f_c(rik) ------------------------------*/
	if(rik < at.R || rik > at.S) drik_fc_rik=0.0;
	if(rik >= at.R && rik <= at.S)
	  drik_fc_rik = -0.5*(PI/(at.S-at.R))*sin(PI*(rik-at.R)/(at.S-at.R));

	cs = (xij*xik + yij*yik + zij*zik)/(rij * rik); /* cos(jik) */
	si = sin(acos(cs));  /* sin(jik) */

	ggg = 1.0 + (c2/d2) - c2/( d2+(at.h-cs)*(at.h-cs) );
	xxx += fc_rik*ggg;

      } /*--- end k loop ---*/

      xxx_n = exp(at.n*log(xxx));  /* xxx^n */
      xxx_n1 = exp((at.n-1.0)*log(xxx));  /* xxx^(n-1) */

      bij = exp( (-1.0/(2.0*at.n)) * log(1.0 + beta_n*xxx_n) );
      bij = at.B * at.x * bij;

      sys.pot += fc_rij*(at.A*exp(-at.lam*rij) - bij*exp(-at.mu*rij)  );

      /*=== force 1 ==========*/
      f1 = drij_fc_rij*(at.A*exp(-at.lam*rij) - bij*exp(-at.mu*rij))
                +fc_rij*((-at.lam)*at.A*exp(-at.lam*rij) 
                                      +(at.mu)*bij*exp(-at.mu*rij));

      sys.fx[i] += -(xij/rij)*f1;
      sys.fy[i] += -(yij/rij)*f1;
      sys.fz[i] += -(zij/rij)*f1;
      sys.fx[j] -= -(xij/rij)*f1;
      sys.fy[j] -= -(yij/rij)*f1;
      sys.fz[j] -= -(zij/rij)*f1;

      /* 2-body part for VIRIAL */
      sys.virX += -(xij/rij)*f1*xij;
      sys.virY += -(yij/rij)*f1*yij;
      sys.virZ += -(zij/rij)*f1*zij;

      for(k=0;k<sys.N;k++) { /*======  2nd k loop =============*/

	if(k==i||k==j) continue;

	xik = sys.rx[i] - sys.rx[k];
	yik = sys.ry[i] - sys.ry[k];
	zik = sys.rz[i] - sys.rz[k];

	/* cyclic boundary condition */
	if(xik > rdx) xik -= sys.Lx;
	if(yik > rdy) yik -= sys.Ly;
	if(zik > rdz) zik -= sys.Lz;
	if(xik < -rdx) xik += sys.Lx;
	if(yik < -rdy) yik += sys.Ly;
	if(zik < -rdz) zik += sys.Lz;

	if(xik > at.S || xik < -at.S) continue;
	if(yik > at.S || yik < -at.S) continue;
	if(zik > at.S || zik < -at.S) continue;

	rik = sqrt(xik*xik + yik*yik + zik*zik);

	/*-----------  f_c(rik) --------------------------------------*/
	fc_rik = 0.0;
	drik_fc_rik = 0.0;

	if(rik > at.S) continue;
	if(rik < at.R) fc_rik = 1.0;
	if(rik >= at.R && rik <= at.S)
	  fc_rik = 0.5 + 0.5 * cos((PI*(rik - at.R))/(at.S - at.R));

	/*-----------  (d/drik)f_c(rik) ------------------------------*/
	if(rik < at.R || rik > at.S) drik_fc_rik=0.0;
	if(rik >= at.R && rik <= at.S)
	  drik_fc_rik = -0.5*(PI/(at.S-at.R))*sin(PI*(rik-at.R)/(at.S-at.R));

	cs = (xij*xik + yij*yik + zij*zik)/(rij * rik); /* cos(jik) */
	si = sin(acos(cs));  /* sin(jik) */

	ggg = 1.0 + (c2/d2) - c2/( d2+(at.h-cs)*(at.h-cs) );

	/*=== force 2 ==========*/
	drik_bij = exp( (-1.0/(2.0*at.n)-1.0) * log(1.0 + beta_n*xxx_n) )
	  * at.B * at.x * (-1.0/(2.0)) * beta_n*xxx_n1*drik_fc_rik*ggg;
	f2 = -fc_rik*drik_bij*exp(-at.mu*rij);

	sys.fx[i] += -(xik/rik)*f2;
	sys.fy[i] += -(yik/rik)*f2;
	sys.fz[i] += -(zik/rik)*f2;

	sys.fx[k] -= -(xik/rik)*f2;
	sys.fy[k] -= -(yik/rik)*f2;
	sys.fz[k] -= -(zik/rik)*f2;

	/* 2-body part for VIRIAL */
	sys.virX += -(xik/rik)*f2*xik;
	sys.virY += -(yik/rik)*f2*yik;
	sys.virZ += -(zik/rik)*f2*zik;

	/*=== force 3 ==========*/
	tmp = d2+(at.h-cs)*(at.h-cs);
	tmp *= tmp;

	dtheta_bij = exp( (-1.0/(2.0*at.n)-1.0) * log(1.0 + beta_n*xxx_n) )
	  * at.B * at.x * (-1.0/(2.0)) * beta_n*xxx_n1
	  * fc_rik*(- (c2*2.0*(at.h-cs)*si)/tmp);

	f3 = -fc_rij*(-dtheta_bij*exp(-at.mu*rij));
	dxj_theta = -(1.0/si)*(1.0/(rij*rik))*(-xik+(rik/rij)*xij*cs);
	dyj_theta = -(1.0/si)*(1.0/(rij*rik))*(-yik+(rik/rij)*yij*cs);
	dzj_theta = -(1.0/si)*(1.0/(rij*rik))*(-zik+(rik/rij)*zij*cs);

	dxk_theta = -(1.0/si)*(1.0/(rij*rik))*(-xij+(rij/rik)*xik*cs);
	dyk_theta = -(1.0/si)*(1.0/(rij*rik))*(-yij+(rij/rik)*yik*cs);
	dzk_theta = -(1.0/si)*(1.0/(rij*rik))*(-zij+(rij/rik)*zik*cs);

	sys.fx[j] += -dxj_theta*f3;
	sys.fy[j] += -dyj_theta*f3;
	sys.fz[j] += -dzj_theta*f3;

	sys.fx[k] += -dxk_theta*f3;
	sys.fy[k] += -dyk_theta*f3;
	sys.fz[k] += -dzk_theta*f3;

	sys.fx[i] -= (-dxj_theta-dxk_theta)*f3;
	sys.fy[i] -= (-dyj_theta-dyk_theta)*f3;
	sys.fz[i] -= (-dzj_theta-dzk_theta)*f3;

	/* 3-body part for VIRIAL */
	sys.virX += (-dxj_theta*xij - dxk_theta*xik)*f3;
	sys.virY += (-dyj_theta*yij - dyk_theta*yik)*f3;
	sys.virZ += (-dzj_theta*zij - dzk_theta*zik)*f3;

      } /*--- end 2nd k loop ---*/

    } /*--- end j loop ---*/
  } /*--- end i loop ---*/

  sys.pot /= 2.0;
  for(i=0;i<sys.N;i++) {
    sys.fx[i] /= 2.0;
    sys.fy[i] /= 2.0;
    sys.fz[i] /= 2.0;
  }

  sys.virX /= 2.0;
  sys.virY /= 2.0;
  sys.virZ /= 2.0;

}
