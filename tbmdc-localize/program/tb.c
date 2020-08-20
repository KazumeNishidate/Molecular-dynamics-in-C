#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/* Tight-Binding Hamiltonian elements */
void   fill_htb1(int local_N, int i, int j, double dr, double s_R, 
		 double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_sssig_r, double v_spsig_r,  
		 double v_ppsig_r, double v_pppi_r)
{
  int m;

  m=0;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+0] = v_sssig_r;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+1] = eru * v_spsig_r;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+2] = emu * v_spsig_r;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+3] = enu * v_spsig_r;

  m=1;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+0] = -eru * v_spsig_r;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+1] = 
    eru2 * v_ppsig_r + (1.0-eru2) * v_pppi_r;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+2] = 
    eru * emu * (v_ppsig_r - v_pppi_r);
  htb.mat1[16*local_N*i+4*j+4*local_N*m+3] = 
    eru * enu * (v_ppsig_r - v_pppi_r);

  m=2;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+0] = -emu * v_spsig_r;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+1] = 
    eru * emu * (v_ppsig_r - v_pppi_r);
  htb.mat1[16*local_N*i+4*j+4*local_N*m+2] = 
    emu2 * v_ppsig_r + (1.0-emu2) * v_pppi_r;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+3] = 
    emu * enu * (v_ppsig_r - v_pppi_r);

  m=3;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+0] = -enu * v_spsig_r;
  htb.mat1[16*local_N*i+4*j+4*local_N*m+1] = 
    eru * enu * (v_ppsig_r - v_pppi_r);
  htb.mat1[16*local_N*i+4*j+4*local_N*m+2] = 
    emu * enu * (v_ppsig_r - v_pppi_r);
  htb.mat1[16*local_N*i+4*j+4*local_N*m+3] = 
    enu2 * v_ppsig_r + (1.0-enu2) * v_pppi_r;

}

/* derivatives of the Tight-Binding Hamiltonian elements */
/* to calculate the Hellman Feymann force                */
void   fill_htb2(int local_N, int i, int j, double dr, double s_R, 
		 double s_dash_R, double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_sssig_r, double v_spsig_r,  
		 double v_ppsig_r, double v_pppi_r, int ii, int jj)
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

  /*
  printf( "distance[%d %d] %f %f %f\n", ii, jj,
          press.dx[ii*sys.N+jj],press.dy[ii*sys.N+jj],press.dz[ii*sys.N+jj] );
  getchar();
  */

/* Fx = d(Htb)/dx */
  m=0;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+0] = ss_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+1] = spx_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+2] = spy_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+3] = spz_dx;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matX[16*local_N*i+4*j+4*local_N*m+0] = ss_dx *press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+1] = spx_dx*press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+2] = spy_dx*press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+3] = spz_dx*press.dx[ii*sys.N+jj];

  m=1;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+0] = pxs_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+1] = pxpx_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+2] = pxpy_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+3] = pxpz_dx;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matX[16*local_N*i+4*j+4*local_N*m+0] = pxs_dx *press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+1] = pxpx_dx*press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+2] = pxpy_dx*press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+3] = pxpz_dx*press.dx[ii*sys.N+jj];

  m=2;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+0] = pys_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+1] = pypx_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+2] = pypy_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+3] = pypz_dx;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matX[16*local_N*i+4*j+4*local_N*m+0] = pys_dx *press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+1] = pypx_dx*press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+2] = pypy_dx*press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+3] = pypz_dx*press.dx[ii*sys.N+jj];

  m=3;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+0] = pzs_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+1] = pzpx_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+2] = pzpy_dx;
  htb.mat2x[16*local_N*i+4*j+4*local_N*m+3] = pzpz_dx;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matX[16*local_N*i+4*j+4*local_N*m+0] = pzs_dx *press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+1] = pzpx_dx*press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+2] = pzpy_dx*press.dx[ii*sys.N+jj];
  press.matX[16*local_N*i+4*j+4*local_N*m+3] = pzpz_dx*press.dx[ii*sys.N+jj];


/* Fy = d(Htb)/dy */
  m=0;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+0] = ss_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+1] = spx_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+2] = spy_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+3] = spz_dy;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matY[16*local_N*i+4*j+4*local_N*m+0] = ss_dy *press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+1] = spx_dy*press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+2] = spy_dy*press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+3] = spz_dy*press.dy[ii*sys.N+jj];

  m=1;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+0] = pxs_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+1] = pxpx_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+2] = pxpy_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+3] = pxpz_dy;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matY[16*local_N*i+4*j+4*local_N*m+0] = pxs_dy *press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+1] = pxpx_dy*press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+2] = pxpy_dy*press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+3] = pxpz_dy*press.dy[ii*sys.N+jj];

  m=2;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+0] = pys_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+1] = pypx_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+2] = pypy_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+3] = pypz_dy;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matY[16*local_N*i+4*j+4*local_N*m+0] = pys_dy *press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+1] = pypx_dy*press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+2] = pypy_dy*press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+3] = pypz_dy*press.dy[ii*sys.N+jj];

  m=3;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+0] = pzs_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+1] = pzpx_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+2] = pzpy_dy;
  htb.mat2y[16*local_N*i+4*j+4*local_N*m+3] = pzpz_dy;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matY[16*local_N*i+4*j+4*local_N*m+0] = pzs_dy *press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+1] = pzpx_dy*press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+2] = pzpy_dy*press.dy[ii*sys.N+jj];
  press.matY[16*local_N*i+4*j+4*local_N*m+3] = pzpz_dy*press.dy[ii*sys.N+jj];


/* Fz = d(Htb)/dz */
  m=0;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+0] = ss_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+1] = spx_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+2] = spy_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+3] = spz_dz;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matZ[16*local_N*i+4*j+4*local_N*m+0] = ss_dz *press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+1] = spx_dz*press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+2] = spy_dz*press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+3] = spz_dz*press.dz[ii*sys.N+jj];

  m=1;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+0] = pxs_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+1] = pxpx_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+2] = pxpy_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+3] = pxpz_dz;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matZ[16*local_N*i+4*j+4*local_N*m+0] = pxs_dz *press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+1] = pxpx_dz*press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+2] = pxpy_dz*press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+3] = pxpz_dz*press.dz[ii*sys.N+jj];

  m=2;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+0] = pys_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+1] = pypx_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+2] = pypy_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+3] = pypz_dz;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matZ[16*local_N*i+4*j+4*local_N*m+0] = pys_dz *press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+1] = pypx_dz*press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+2] = pypy_dz*press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+3] = pypz_dz*press.dz[ii*sys.N+jj];

  m=3;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+0] = pzs_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+1] = pzpx_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+2] = pzpy_dz;
  htb.mat2z[16*local_N*i+4*j+4*local_N*m+3] = pzpz_dz;

  /* --- taka (Wednesday, June 7, 2000) --- */
  press.matZ[16*local_N*i+4*j+4*local_N*m+0] = pzs_dz *press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+1] = pzpx_dz*press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+2] = pzpy_dz*press.dz[ii*sys.N+jj];
  press.matZ[16*local_N*i+4*j+4*local_N*m+3] = pzpz_dz*press.dz[ii*sys.N+jj];

}

void   diagonal(int local_N, int i) {

  int ii, jj;

  htb.mat1[16*local_N*i+4*i+4*local_N*0+0] = htb.e_s;
  htb.mat1[16*local_N*i+4*i+4*local_N*0+1] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*0+2] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*0+3] = 0.0;

  htb.mat1[16*local_N*i+4*i+4*local_N*1+0] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*1+1] = htb.e_p;
  htb.mat1[16*local_N*i+4*i+4*local_N*1+2] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*1+3] = 0.0;

  htb.mat1[16*local_N*i+4*i+4*local_N*2+0] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*2+1] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*2+2] = htb.e_p;
  htb.mat1[16*local_N*i+4*i+4*local_N*2+3] = 0.0;

  htb.mat1[16*local_N*i+4*i+4*local_N*3+0] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*3+1] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*3+2] = 0.0;
  htb.mat1[16*local_N*i+4*i+4*local_N*3+3] = htb.e_p;

  for(ii=0;ii<4;ii++) {
    for(jj=0;jj<4;jj++) {
      htb.mat2x[16*local_N*i+4*i+4*local_N*ii+jj] = 0.0;
      htb.mat2y[16*local_N*i+4*i+4*local_N*ii+jj] = 0.0;
      htb.mat2z[16*local_N*i+4*i+4*local_N*ii+jj] = 0.0;
    }
  }
}




