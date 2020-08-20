#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "headers/md.h"

void   set_Vps_table(int n_pl, double *kx, double *ky, double *kz)
{
  int nq1, nq2, type;

  /*--- BHS psuedo potential calculation ---*/
  for(type=0;type<atom.ntypes;type++) {
    for(nq1=0;nq1<n_pl;nq1++) {
      for(nq2=0;nq2<n_pl;nq2++) {

	calc_core_Vps(type, nq1, nq2, kx, ky, kz);
	
	if(nq1==0 || nq2 ==0) {
	  calc_delta_Vps_zero(type, nq1, nq2, kx, ky, kz);
	} else {
	  calc_delta_Vps(type, nq1, nq2, kx, ky, kz);
	}

	wk.Vps /= cell.vol;
	wk.Vpsmat[type][nq1][nq2] = wk.Vps;
      }
    }
  }

}


void   set_Vps(int n_pl, double *kx, double *ky, double *kz)
{
  int i, nq1, nq2, ib, type;
  double dgx, dgy, dgz, ggr, psr, psi;

  /*--- set up H*C with Vps ---*/  
  for(nq1=0;nq1<n_pl;nq1++) {

    for(ib=0;ib<atom.nvband_all;ib++) { /* initialize */
      wv.hcgr[ib][nq1] = 0.0;
      wv.hcgi[ib][nq1] = 0.0;
    }
    for(nq2=0;nq2<n_pl;nq2++) {
      dgx = kx[nq2]-kx[nq1];
      dgy = ky[nq2]-ky[nq1];
      dgz = kz[nq2]-kz[nq1];
      psr = 0.0;
      psi = 0.0;
      for(i=0;i<atom.N;i++) {
	type = atom.type[i];
	ggr = dgx*cell.rx[i]+dgy*cell.ry[i]+dgz*cell.rz[i];
	psr += cos(ggr)*(wk.Vpsmat[type][nq1][nq2]);
	psi += sin(ggr)*(wk.Vpsmat[type][nq1][nq2]);
      }
      for(ib=0;ib<atom.nvband_all;ib++) {
	wv.hcgr[ib][nq1] += psr*wv.cg2r[ib][nq2] - psi*wv.cg2i[ib][nq2];
	wv.hcgi[ib][nq1] += psi*wv.cg2r[ib][nq2] + psr*wv.cg2i[ib][nq2];
      }
    }

  }
  
}

void  clear_hmat_small()
{
  int nq1, nq2, cnt=0;
  
  for(nq1=0;nq1<wv.nplwin;nq1++)
    for(nq2=0;nq2<=nq1;nq2++) {
      hm.mat_re[cnt] = 0.0;
      hm.mat_im[cnt] = 0.0;
      cnt++;
    }
}

void   set_Ek(int npw, double *kx, double *ky, double *kz)
{
  int nq1, nq2, cnt=0;
  double gx, gy, gz;

  for(nq1=0;nq1<npw;nq1++) {
    for(nq2=0;nq2<=nq1;nq2++) {
      
      if(nq1==nq2) { /* electronic kinetic energy (1/2) |k+G|^2 */
	gx = kx[nq1];
	gy = ky[nq1];
	gz = kz[nq1];
	hm.mat_re[cnt] += 0.5*(gx*gx+gy*gy+gz*gz);
      }
      cnt++;
    }
  }
}

void   set_Vps_hmat(int n_pl, double *kx, double *ky, double *kz)
{
  int i, nq1, nq2, type, cnt;
  double dgx, dgy, dgz, ggr;
  
  /*--- BHS psuedo potential calculation ---*/
  for(type=0;type<atom.ntypes;type++) {
    for(nq1=0;nq1<n_pl;nq1++) {
      for(nq2=0;nq2<=nq1;nq2++) {
	
	calc_core_Vps(type, nq1, nq2, kx, ky, kz);
	
	if(nq1==0 || nq2 ==0) {
	  calc_delta_Vps_zero(type, nq1, nq2, kx, ky, kz);
	} else {
	  calc_delta_Vps(type, nq1, nq2, kx, ky, kz);
	}

	wk.Vps /= cell.vol;
	wk.Vpsmat[type][nq1][nq2] = wk.Vps;
      }
    }
  }
  
  /*--- fill up the Hamiltonian with Vps to diagonalize ---*/
  cnt=0;
  for(nq1=0;nq1<n_pl;nq1++) {
    for(nq2=0;nq2<=nq1;nq2++) {

      dgx = kx[nq2]-kx[nq1];
      dgy = ky[nq2]-ky[nq1];
      dgz = kz[nq2]-kz[nq1];

      for(i=0;i<atom.N;i++) { 
	type = atom.type[i];
	ggr = dgx*cell.rx[i]+dgy*cell.ry[i]+dgz*cell.rz[i];
	hm.mat_re[cnt] += cos(ggr)*(wk.Vpsmat[type][nq1][nq2]);
	hm.mat_im[cnt] += sin(ggr)*(wk.Vpsmat[type][nq1][nq2]);
      }
      cnt++;
    }
  }
  
}

void   mk_Bessel(int n_pl, double *kx, double *ky, double *kz)
{                              /* making up the Bessel function table */
  int eru, nq1, nq2;
  int type;
  
  double k1, k2, k1d, k2d;
  double z1, z2, z3, fz1, fz2, fz3;
  double bm1, bm2, bm3, dbm1, dbm2, dbm3;
  double e1, e2, e3;
  double s1, s2, s3;
  double c11, c12, c13, c21, c22, c23, c31, c32, c33;
  double ds1, ds2, ds3;
  double alpha[4], aa[7];
  
  for(type=0;type<atom.ntypes;type++) {
    
    for(nq1=1;nq1<n_pl;nq1++) {
      k1d = kx[nq1]*kx[nq1]+ky[nq1]*ky[nq1]+kz[nq1]*kz[nq1];
      k1  = sqrt(k1d);

      for(nq2=1;nq2<=nq1;nq2++) { 
	k2d = kx[nq2]*kx[nq2]+ky[nq2]*ky[nq2]+kz[nq2]*kz[nq2];
	k2  = sqrt(k2d);
	
	for(eru=0;eru<atom.norbit[type];eru++) {
	  alpha[1] = bhs.alpha[type][eru][1];
	  alpha[2] = bhs.alpha[type][eru][2];
	  alpha[3] = bhs.alpha[type][eru][3];
	  aa[1] = bhs.A[type][eru][1];
	  aa[2] = bhs.A[type][eru][2];
	  aa[3] = bhs.A[type][eru][3];
	  aa[4] = bhs.A[type][eru][4];
	  aa[5] = bhs.A[type][eru][5];
	  aa[6] = bhs.A[type][eru][6];
	  
	  z1 = (k1*k2)/(2.0*alpha[1]);
	  z2 = (k1*k2)/(2.0*alpha[2]);
	  z3 = (k1*k2)/(2.0*alpha[3]);
	  fz1 = sqrt(2.0/(PI*z1));
	  fz2 = sqrt(2.0/(PI*z2));
	  fz3 = sqrt(2.0/(PI*z3));
	  
	  switch (eru) {
	  case 0:
	    bm1  = fz1*sinh(z1);
	    bm2  = fz2*sinh(z2);
	    bm3  = fz3*sinh(z3);
	    dbm1 = fz1*( cosh(z1)-sinh(z1)/(2.0*z1) );
	    dbm2 = fz2*( cosh(z2)-sinh(z2)/(2.0*z2) );
	    dbm3 = fz3*( cosh(z3)-sinh(z3)/(2.0*z3) );
	    break;
	  case 1:
	    bm1  = fz1*( cosh(z1)-sinh(z1)/z1 );
	    bm2  = fz2*( cosh(z2)-sinh(z2)/z2 );
	    bm3  = fz3*( cosh(z3)-sinh(z3)/z3 );
	    dbm1 = fz1*sinh(z1) - bm1*3.0/(2.0*z1);
	    dbm2 = fz2*sinh(z2) - bm2*3.0/(2.0*z2);
	    dbm3 = fz3*sinh(z3) - bm3*3.0/(2.0*z3);
	    break;
	  case 2:
	    bm1  = fz1*( (1.0+3.0/(z1*z1))*sinh(z1) - 3.0*cosh(z1)/z1 );
	    bm2  = fz2*( (1.0+3.0/(z2*z2))*sinh(z2) - 3.0*cosh(z2)/z2 );
	    bm3  = fz3*( (1.0+3.0/(z3*z3))*sinh(z3) - 3.0*cosh(z3)/z3 );
	    dbm1 = fz1*( cosh(z1)-sinh(z1)/z1 ) - bm1*5.0/(2.0*z1);
	    dbm2 = fz2*( cosh(z2)-sinh(z2)/z2 ) - bm2*5.0/(2.0*z2);
	    dbm3 = fz3*( cosh(z3)-sinh(z3)/z3 ) - bm3*5.0/(2.0*z3);
	    break;
	  case 3:
	    bm1  = fz1*( (1.0+15.0/(z1*z1))*cosh(z1) -
			 (6.0/z1+15.0/(z1*z1*z1))*sinh(z1) );
	    bm2  = fz2*( (1.0+15.0/(z2*z2))*cosh(z2) -
			 (6.0/z2+15.0/(z2*z2*z2))*sinh(z2) );
	    bm3  = fz3*( (1.0+15.0/(z3*z3))*cosh(z3) -
			 (6.0/z3+15.0/(z3*z3*z3))*sinh(z3) );
	    dbm1 = fz1*( (1.0+3.0/(z1*z1))*sinh(z1) - 3.0*cosh(z1)/z1 ) -
	      bm1*7.0/(2.0*z1);
	    dbm2 = fz2*( (1.0+3.0/(z2*z2))*sinh(z2) - 3.0*cosh(z2)/z2 ) -
	      bm2*7.0/(2.0*z2);
	    dbm3 = fz3*( (1.0+3.0/(z3*z3))*sinh(z3) - 3.0*cosh(z3)/z3 ) -
	      bm3*7.0/(2.0*z3);
	    break;
	  default:
	    bm1=bm2=bm3=dbm1=dbm2=dbm3=0.0;
	    printf("(mkbessel) eru = %d is not supported\n",eru);
	    break;
	  }
	  
	  e1 = exp( -(k1d+k2d)/(4.0*alpha[1]) );
	  e2 = exp( -(k1d+k2d)/(4.0*alpha[2]) );
	  e3 = exp( -(k1d+k2d)/(4.0*alpha[3]) );
	  
	  s1 = e1*bm1*aa[1]/alpha[1];
	  s2 = e2*bm2*aa[2]/alpha[2];
	  s3 = e3*bm3*aa[3]/alpha[3];
	  
	  c11 = -1.0/alpha[1];
	  c21 = -1.0/alpha[2];
	  c31 = -1.0/alpha[3];
	  
	  c12 = (k1d+k2d)/(4.0*alpha[1]*alpha[1]);
	  c22 = (k1d+k2d)/(4.0*alpha[2]*alpha[2]);
	  c32 = (k1d+k2d)/(4.0*alpha[3]*alpha[3]);
	  
	  c13 = -(k1*k2)/(2.0*alpha[1]*alpha[1]);
	  c23 = -(k1*k2)/(2.0*alpha[2]*alpha[2]);
	  c33 = -(k1*k2)/(2.0*alpha[3]*alpha[3]);
	  
	  ds1 = -(aa[4]/alpha[1])*(c11*e1*bm1+c12*e1*bm1+c13*e1*dbm1);
	  ds2 = -(aa[5]/alpha[2])*(c21*e2*bm2+c22*e2*bm2+c23*e2*dbm2);
	  ds3 = -(aa[6]/alpha[3])*(c31*e3*bm3+c32*e3*bm3+c33*e3*dbm3);
	  
	  bhs.bsl[eru][type][bhs.id[nq1]][bhs.id[nq2]] =
	    ( PI/(4.0*sqrt(k1*k2)) )*( s1+s2+s3+ds1+ds2+ds3 );
	  bhs.bsl[eru][type][bhs.id[nq2]][bhs.id[nq1]] =
	    bhs.bsl[eru][type][bhs.id[nq1]][bhs.id[nq2]];
	}
	
      }
    }
    
  }

}

void calc_delta_Vps(int type, int nq1, int nq2,
		    double *kx, double *ky, double *kz)
{
  double kx1, ky1, kz1, k1;
  double kx2, ky2, kz2, k2;
  double angle, angle2, angle3, P0, P1, P2, P3;

  kx1 = kx[nq1];  ky1 = ky[nq1];  kz1 = kz[nq1];
  k1  = sqrt(kx1*kx1+ky1*ky1+kz1*kz1);
  kx2 = kx[nq2];  ky2 = ky[nq2];  kz2 = kz[nq2];
  k2  = sqrt(kx2*kx2+ky2*ky2+kz2*kz2);

  angle = (kx1*kx2+ky1*ky2+kz1*kz2)/(k1*k2);
  angle2 = angle*angle;
  angle3 = angle*angle*angle;
  P0 = 1.0;
  P1 = angle;
  P2 = 0.5*(3.0*angle2-1.0);
  P3 = 0.5*(5.0*angle3-3.0*angle);

  wk.Vps += 4.0*PI*( (1.0)*P0*bhs.bsl[0][type][bhs.id[nq1]][bhs.id[nq2]]+
                 (2.0+1.0)*P1*bhs.bsl[1][type][bhs.id[nq1]][bhs.id[nq2]]+
                 (4.0+1.0)*P2*bhs.bsl[2][type][bhs.id[nq1]][bhs.id[nq2]]+
		 (6.0+1.0)*P3*bhs.bsl[3][type][bhs.id[nq1]][bhs.id[nq2]] );

}

void calc_delta_Vps_zero(int type, int nq1, int nq2,
			 double *kx, double *ky, double *kz)
{
  double kx1, ky1, kz1, k1d;
  double kx2, ky2, kz2, k2d;
  double a1, a2, a3, z;
  double c1, c2, c3, c4, c5, c6;
  double a11, a21, a31, a12, a22, a32;

  kx1 = kx[nq1];  ky1 = ky[nq1];  kz1 = kz[nq1];
  k1d = kx1*kx1+ky1*ky1+kz1*kz1;

  kx2 = kx[nq2];  ky2 = ky[nq2];  kz2 = kz[nq2];
  k2d = kx2*kx2+ky2*ky2+kz2*kz2;

  a1 = bhs.alpha[type][0][1];
  a2 = bhs.alpha[type][0][2];
  a3 = bhs.alpha[type][0][3];

  /* x^y = exp(y*log(x))  [y = integer, x><0] or [y=float, x>0]  */
  a11 = exp(1.5*log(a1));  a21 = exp(1.5*log(a2));
  a31 = exp(1.5*log(a3));  a12 = exp(2.5*log(a1));
  a22 = exp(2.5*log(a2));  a32 = exp(2.5*log(a3));

  c1 = bhs.A[type][0][1];  c2 = bhs.A[type][0][2];
  c3 = bhs.A[type][0][3];  c4 = bhs.A[type][0][4];
  c5 = bhs.A[type][0][5];  c6 = bhs.A[type][0][6];

  z = -(k1d+k2d)/4.0;

  wk.Vps += PI*SQRTPI*( 
               (c1/a11)*exp(z/a1) 
             + (c2/a21)*exp(z/a2)
             + (c3/a31)*exp(z/a3) 
             + (c4/a12)*(1.5+z/a1)*exp(z/a1)
             + (c5/a22)*(1.5+z/a2)*exp(z/a2)
             + (c6/a32)*(1.5+z/a3)*exp(z/a3) );

}

void calc_core_Vps(int type, int nq1, int nq2, 
		   double *kx, double *ky, double *kz)
{
  double a1, a2, c1, c2;
  double dkx, dky, dkz, gnorm2;

  c1 = bhs.core_orbit_c[type][0];
  c2 = bhs.core_orbit_c[type][1];
  a1 = bhs.core_orbit_alpha[type][0];
  a2 = bhs.core_orbit_alpha[type][1];

  if(nq1==nq2) {
    wk.Vps = PI*atom.Zv[type]*( c1/a1 + c2/a2 );  /* BUG (natsu) */
  } else {

    dkx = kx[nq1]-kx[nq2];
    dky = ky[nq1]-ky[nq2];
    dkz = kz[nq1]-kz[nq2];
    gnorm2 = dkx*dkx+dky*dky+dkz*dkz;

    wk.Vps = -(4.0*PI*atom.Zv[type]/gnorm2)*
      ( c1*exp(-gnorm2/(4.0*a1)) + c2*exp(-gnorm2/(4.0*a2)) );

  }

}

void   get_HF_force(int n_pl, double *kx, double *ky, double *kz)
{
  int i, ib, nq1, nq2, type;
  double ggr;
  double dgx, dgy, dgz;
  double vc, vs, fi;

  for(i=0;i<atom.N;i++) {
    wk.fecx[i] = 0.0;
    wk.fecy[i] = 0.0;
    wk.fecz[i] = 0.0;
  }

  for(i=0;i<atom.N;i++) {        
    type = atom.type[i];

    for(nq1=0;nq1<n_pl;nq1++) {
      for(nq2=0;nq2<nq1;nq2++) {  /* only left-down part */

	dgx = kx[nq2]-kx[nq1];
	dgy = ky[nq2]-ky[nq1];
	dgz = kz[nq2]-kz[nq1];
	ggr = dgx*cell.rx[i]+dgy*cell.ry[i]+dgz*cell.rz[i];
	vc = cos(ggr);
	vs = sin(ggr);

	for(ib=0;ib<atom.nvband_all;ib++) {/*-------------------------*/
	  fi = ( (wv.cg2r[ib][nq1]*wv.cg2i[ib][nq2] 
		- wv.cg2i[ib][nq1]*wv.cg2r[ib][nq2])*vc
	       + (wv.cg2r[ib][nq1]*wv.cg2r[ib][nq2]
	        + wv.cg2i[ib][nq1]*wv.cg2i[ib][nq2])*vs )
	    *wk.Vpsmat[type][nq1][nq2]*2.0*atom.fermi[ib];

	  wk.fecx[i] += dgx*fi*2.0;
	  wk.fecy[i] += dgy*fi*2.0;
	  wk.fecz[i] += dgz*fi*2.0;
	} /*----------------------------------------------------------*/

      }
    }

  }

}
