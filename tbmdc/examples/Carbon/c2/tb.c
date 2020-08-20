#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/* Tight-Binding Hamiltonian elements */
void   fill_htb1(int i, int j, double dr, double s_R, 
		 double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_sssig_r, double v_spsig_r,  
		 double v_ppsig_r, double v_pppi_r)
{
  int m;
  int NN; /* check? : taka ( May 31, 2000 ) */
  m=0;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat1[NN+0] =  v_sssig_r;
  htb.mat1[NN+1] =  eru * v_spsig_r;
  htb.mat1[NN+2] =  emu * v_spsig_r;
  htb.mat1[NN+3] =  enu * v_spsig_r;

  m=1;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat1[NN+0] = -eru * v_spsig_r;
  htb.mat1[NN+1] =  eru2 * v_ppsig_r + (1.0-eru2) * v_pppi_r;
  htb.mat1[NN+2] =  eru * emu * (v_ppsig_r - v_pppi_r);
  htb.mat1[NN+3] =  eru * enu * (v_ppsig_r - v_pppi_r);

  m=2;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat1[NN+0] = -emu * v_spsig_r;
  htb.mat1[NN+1] =  eru * emu * (v_ppsig_r - v_pppi_r);
  htb.mat1[NN+2] =  emu2 * v_ppsig_r + (1.0-emu2) * v_pppi_r;
  htb.mat1[NN+3] =  emu * enu * (v_ppsig_r - v_pppi_r);

  m=3;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat1[NN+0] = -enu * v_spsig_r;
  htb.mat1[NN+1] =  eru * enu * (v_ppsig_r - v_pppi_r);
  htb.mat1[NN+2] =  emu * enu * (v_ppsig_r - v_pppi_r);
  htb.mat1[NN+3] =  enu2 * v_ppsig_r + (1.0-enu2) * v_pppi_r;
}

/* derivatives of the Tight-Binding Hamiltonian elements */
/* to calculate the Hellman Feymann force                */
void   fill_htb2(int i, int j, double dr, double s_R, double s_dash_R,
	 double eru, double emu, double enu,
	 double eru2, double emu2, double enu2, 
	 double v_sssig_r, double v_spsig_r,  
	 double v_ppsig_r, double v_pppi_r)
{
  int m;
  double ss_dx, spx_dx, spy_dx, spz_dx;
  double pxs_dx, pxpx_dx, pxpy_dx, pxpz_dx;
  double pys_dx, pypx_dx, pypy_dx, pypz_dx;
  double pzs_dx, pzpx_dx, pzpy_dx, pzpz_dx;
  double ss_dy, spx_dy, spy_dy, spz_dy;
  double pxs_dy, pxpx_dy, pxpy_dy, pxpz_dy;
  double pys_dy, pypx_dy, pypy_dy, pypz_dy;
  double pzs_dy, pzpx_dy, pzpy_dy, pzpz_dy;
  double ss_dz, spx_dz, spy_dz, spz_dz;
  double pxs_dz, pxpx_dz, pxpy_dz, pxpz_dz;
  double pys_dz, pypx_dz, pypy_dz, pypz_dz;
  double pzs_dz, pzpx_dz, pzpy_dz, pzpz_dz;

  int NN; /* check? : taka ( May 31, 2000 ) */

/* Fx = [d(T_ab)/dx] */
  ss_dx  = eru * htb.v_sssig * s_dash_R;
  spx_dx = (1.0/dr)*((1.0-eru2)*v_spsig_r + eru2*dr*htb.v_spsig*s_dash_R);
  spy_dx = (eru*emu/dr)*(-v_spsig_r + dr*htb.v_spsig*s_dash_R);
  spz_dx = (eru*enu/dr)*(-v_spsig_r + dr*htb.v_spsig*s_dash_R);

  pxs_dx  = -spx_dx;
  pxpx_dx = (eru/dr)*(2.0*(1.0-eru2)*(v_ppsig_r - v_pppi_r)
   +eru2*dr*htb.v_ppsig*s_dash_R+(1.0-eru2)*dr*htb.v_pppi*s_dash_R);
  pxpy_dx = (emu/dr)*((1.0-2.0*eru2)*(v_ppsig_r - v_pppi_r)
   +eru2*dr*(htb.v_ppsig*s_dash_R-htb.v_pppi*s_dash_R));
  pxpz_dx = (enu/dr)*((1.0-2.0*eru2)*(v_ppsig_r - v_pppi_r)
   +eru2*dr*(htb.v_ppsig*s_dash_R-htb.v_pppi*s_dash_R));

  pys_dx  = -spy_dx;
  pypx_dx = pxpy_dx;
  pypy_dx = (eru/dr)*(-2.0*emu2*(v_ppsig_r - v_pppi_r)
   +emu2*dr*htb.v_ppsig*s_dash_R+(1.0-emu2)*dr*htb.v_pppi*s_dash_R);
  pypz_dx = (eru*emu*enu/dr)*(-2.0*(v_ppsig_r - v_pppi_r)
   +dr*(htb.v_ppsig*s_dash_R-htb.v_pppi*s_dash_R));

  pzs_dx  = -spz_dx;
  pzpx_dx = pxpz_dx;
  pzpy_dx = pypz_dx;
  pzpz_dx = (eru/dr)*(-2.0*enu2*(v_ppsig_r - v_pppi_r)
   +enu2*dr*htb.v_ppsig*s_dash_R+(1.0-enu2)*dr*htb.v_pppi*s_dash_R);

/* Fy = [d(T_ab)/dy] */
  ss_dy  = emu * htb.v_sssig * s_dash_R;
  spx_dy = spy_dx;
  spy_dy = (1.0/dr)*((1.0-emu2)*v_spsig_r + emu2*dr*htb.v_spsig*s_dash_R);
  spz_dy = (emu*enu/dr)*(-v_spsig_r + dr*htb.v_spsig*s_dash_R);

  pxs_dy  = -spx_dy;
  pxpx_dy = (emu/dr)*(-2.0*eru2*(v_ppsig_r - v_pppi_r)
   +eru2*dr*htb.v_ppsig*s_dash_R+(1.0-eru2)*dr*htb.v_pppi*s_dash_R);
  pxpy_dy = (eru/dr)*((1.0-2.0*emu2)*(v_ppsig_r - v_pppi_r)
   +emu2*dr*(htb.v_ppsig*s_dash_R-htb.v_pppi*s_dash_R));
  pxpz_dy = pypz_dx;

  pys_dy  = -spy_dy;
  pypx_dy = pxpy_dy;
  pypy_dy = (emu/dr)*(2.0*(1.0-emu2)*(v_ppsig_r - v_pppi_r)
   +emu2*dr*htb.v_ppsig*s_dash_R+(1.0-emu2)*dr*htb.v_pppi*s_dash_R);
  pypz_dy = (enu/dr)*((1.0-2.0*emu2)*(v_ppsig_r - v_pppi_r)
   +emu2*dr*(htb.v_ppsig*s_dash_R-htb.v_pppi*s_dash_R));

  pzs_dy  = -spz_dy;
  pzpx_dy = pxpz_dy;
  pzpy_dy = pypz_dy;
  pzpz_dy = (emu/dr)*(-2.0*enu2*(v_ppsig_r - v_pppi_r)
   +enu2*dr*htb.v_ppsig*s_dash_R+(1.0-enu2)*dr*htb.v_pppi*s_dash_R);

/* Fz = [d(T_ab)/dz] */
  ss_dz  = enu * htb.v_sssig * s_dash_R;
  spx_dz = (eru*enu/dr)*(-v_spsig_r + dr*htb.v_spsig*s_dash_R);
  spy_dz = (emu*enu/dr)*(-v_spsig_r + dr*htb.v_spsig*s_dash_R);
  spz_dz = (1.0/dr)*((1.0-enu2)*v_spsig_r + enu2*dr*htb.v_spsig*s_dash_R);

  pxs_dz  = -spx_dz;
  pxpx_dz = (enu/dr)*(-2.0*eru2*(v_ppsig_r - v_pppi_r)
   +eru2*dr*htb.v_ppsig*s_dash_R+(1.0-eru2)*dr*htb.v_pppi*s_dash_R);
  pxpy_dz = pypz_dx;
  pxpz_dz = (eru/dr)*((1.0-2.0*enu2)*(v_ppsig_r - v_pppi_r)
   +enu2*dr*(htb.v_ppsig*s_dash_R-htb.v_pppi*s_dash_R));

  pys_dz  = -spy_dz;
  pypx_dz = pxpy_dz;
  pypy_dz = (enu/dr)*(-2.0*emu2*(v_ppsig_r - v_pppi_r)
   +emu2*dr*htb.v_ppsig*s_dash_R+(1.0-emu2)*dr*htb.v_pppi*s_dash_R);
  pypz_dz = (emu/dr)*((1.0-2.0*enu2)*(v_ppsig_r - v_pppi_r)
   +enu2*dr*(htb.v_ppsig*s_dash_R-htb.v_pppi*s_dash_R));

  pzs_dz  = -spz_dz;
  pzpx_dz = pxpz_dz;
  pzpy_dz = pypz_dz;
  pzpz_dz = (enu/dr)*(2.0*(1.0-enu2)*(v_ppsig_r - v_pppi_r)
   +enu2*dr*htb.v_ppsig*s_dash_R+(1.0-enu2)*dr*htb.v_pppi*s_dash_R);

/* Fx = d(Htb)/dx */
  m=0;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2x[NN+0] = ss_dx;
  htb.mat2x[NN+1] = spx_dx;
  htb.mat2x[NN+2] = spy_dx;
  htb.mat2x[NN+3] = spz_dx;

  m=1;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2x[NN+0] = pxs_dx;
  htb.mat2x[NN+1] = pxpx_dx;
  htb.mat2x[NN+2] = pxpy_dx;
  htb.mat2x[NN+3] = pxpz_dx;

  m=2;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2x[NN+0] = pys_dx;
  htb.mat2x[NN+1] = pypx_dx;
  htb.mat2x[NN+2] = pypy_dx;
  htb.mat2x[NN+3] = pypz_dx;

  m=3;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2x[NN+0] = pzs_dx;
  htb.mat2x[NN+1] = pzpx_dx;
  htb.mat2x[NN+2] = pzpy_dx;
  htb.mat2x[NN+3] = pzpz_dx;

/* Fy = d(Htb)/dy */
  m=0;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2y[NN+0] = ss_dy;
  htb.mat2y[NN+1] = spx_dy;
  htb.mat2y[NN+2] = spy_dy;
  htb.mat2y[NN+3] = spz_dy;

  m=1;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2y[NN+0] = pxs_dy;
  htb.mat2y[NN+1] = pxpx_dy;
  htb.mat2y[NN+2] = pxpy_dy;
  htb.mat2y[NN+3] = pxpz_dy;

  m=2;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2y[NN+0] = pys_dy;
  htb.mat2y[NN+1] = pypx_dy;
  htb.mat2y[NN+2] = pypy_dy;
  htb.mat2y[NN+3] = pypz_dy;

  m=3;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2y[NN+0] = pzs_dy;
  htb.mat2y[NN+1] = pzpx_dy;
  htb.mat2y[NN+2] = pzpy_dy;
  htb.mat2y[NN+3] = pzpz_dy;

/* Fz = d(Htb)/dz */
  m=0;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2z[NN+0] = ss_dz;
  htb.mat2z[NN+1] = spx_dz;
  htb.mat2z[NN+2] = spy_dz;
  htb.mat2z[NN+3] = spz_dz;

  m=1;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2z[NN+0] = pxs_dz;
  htb.mat2z[NN+1] = pxpx_dz;
  htb.mat2z[NN+2] = pxpy_dz;
  htb.mat2z[NN+3] = pxpz_dz;

  m=2;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2z[NN+0] = pys_dz;
  htb.mat2z[NN+1] = pypx_dz;
  htb.mat2z[NN+2] = pypy_dz;
  htb.mat2z[NN+3] = pypz_dz;

  m=3;
  NN = 16*sys.N*i+4*j+4*sys.N*m;
  htb.mat2z[NN+0] = pzs_dz;
  htb.mat2z[NN+1] = pzpx_dz;
  htb.mat2z[NN+2] = pzpy_dz;
  htb.mat2z[NN+3] = pzpz_dz;
}

void   diagonal(int i) {

  int ii, jj;

  int NN; /* check? : taka ( May 31, 2000 ) */

  NN = 16*sys.N*i+4*i+4*sys.N*0;
  htb.mat1[NN+0] = htb.e_s;  htb.mat1[NN+1] = 0.0;
  htb.mat1[NN+2] = 0.0;      htb.mat1[NN+3] = 0.0;

  NN = 16*sys.N*i+4*i+4*sys.N*1;
  htb.mat1[NN+0] = 0.0;      htb.mat1[NN+1] = htb.e_p;
  htb.mat1[NN+2] = 0.0;      htb.mat1[NN+3] = 0.0;

  NN = 16*sys.N*i+4*i+4*sys.N*2;
  htb.mat1[NN+0] = 0.0;      htb.mat1[NN+1] = 0.0;
  htb.mat1[NN+2] = htb.e_p;  htb.mat1[NN+3] = 0.0;

  NN = 16*sys.N*i+4*i+4*sys.N*3;
  htb.mat1[NN+0] = 0.0;      htb.mat1[NN+1] = 0.0;
  htb.mat1[NN+2] = 0.0;      htb.mat1[NN+3] = htb.e_p;

  for(ii=0;ii<4;ii++) {
    for(jj=0;jj<4;jj++) {
      NN = 16*sys.N*i+4*i+4*sys.N*ii+jj;
      htb.mat2x[NN] = 0.0;  htb.mat2y[NN] = 0.0;  htb.mat2z[NN] = 0.0;
    }
  }
}

