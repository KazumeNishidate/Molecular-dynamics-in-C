#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*--------------------------------------------------------------------------*/
/* Tight-Binding Hamiltonian elements                                       */
/*--------------------------------------------------------------------------*/
void   fill_htb1(int i, int j, double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_A_r, double v_B_r, double v_C_r, double v_D_r)
{
  htb.mat1[16*sys.N*i+4*j+4*sys.N*0+0] = v_A_r;
  htb.mat1[16*sys.N*i+4*j+4*sys.N*0+1] = eru*v_B_r;
  htb.mat1[16*sys.N*i+4*j+4*sys.N*0+2] = emu*v_B_r;
  htb.mat1[16*sys.N*i+4*j+4*sys.N*0+3] = enu*v_B_r;

  htb.mat1[16*sys.N*i+4*j+4*sys.N*1+0] = -eru*v_B_r;
  htb.mat1[16*sys.N*i+4*j+4*sys.N*1+1] = eru2*v_C_r + (1.0 - eru2)*v_D_r;
  htb.mat1[16*sys.N*i+4*j+4*sys.N*1+2] = eru*emu*(v_C_r - v_D_r);
  htb.mat1[16*sys.N*i+4*j+4*sys.N*1+3] = eru*enu*(v_C_r - v_D_r);

  htb.mat1[16*sys.N*i+4*j+4*sys.N*2+0] = -emu*v_B_r;
  htb.mat1[16*sys.N*i+4*j+4*sys.N*2+1] = eru*emu*(v_C_r - v_D_r);
  htb.mat1[16*sys.N*i+4*j+4*sys.N*2+2] = emu2*v_C_r + (1.0 - emu2)*v_D_r;
  htb.mat1[16*sys.N*i+4*j+4*sys.N*2+3] = emu*enu*(v_C_r - v_D_r);

  htb.mat1[16*sys.N*i+4*j+4*sys.N*3+0] = -enu*v_B_r;
  htb.mat1[16*sys.N*i+4*j+4*sys.N*3+1] = eru*enu*(v_C_r - v_D_r);
  htb.mat1[16*sys.N*i+4*j+4*sys.N*3+2] = emu*enu*(v_C_r - v_D_r);
  htb.mat1[16*sys.N*i+4*j+4*sys.N*3+3] = enu2*v_C_r + (1.0 - enu2)*v_D_r;
}

/*--------------------------------------------------------------------------*/
/* derivatives of the Tight-Binding Hamiltonian elements to calculate       */
/* the Hellman Feymann force                                                */
/*--------------------------------------------------------------------------*/
void   fill_htb2(int i, int j, double dr,
		 double ds_R_A, double ds_R_B, double ds_R_C, double ds_R_D,
		 double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_A_r, double v_B_r, double v_C_r, double v_D_r)
{
  double  ss_dx,  spx_dx,  spy_dx,  spz_dx;
  double pxs_dx, pxpx_dx, pxpy_dx, pxpz_dx;
  double pys_dx, pypx_dx, pypy_dx, pypz_dx;
  double pzs_dx, pzpx_dx, pzpy_dx, pzpz_dx;

  double  ss_dy,  spx_dy,  spy_dy,  spz_dy;
  double pxs_dy, pxpx_dy, pxpy_dy, pxpz_dy;
  double pys_dy, pypx_dy, pypy_dy, pypz_dy;
  double pzs_dy, pzpx_dy, pzpy_dy, pzpz_dy;

  double  ss_dz,  spx_dz,  spy_dz,  spz_dz;
  double pxs_dz, pxpx_dz, pxpy_dz, pxpz_dz;
  double pys_dz, pypx_dz, pypy_dz, pypz_dz;
  double pzs_dz, pzpx_dz, pzpy_dz, pzpz_dz;

/* Fx = [d(T_ab)/dx] */
  ss_dx  = eru*htb.v_A*ds_R_A;
  spx_dx = (1.0/dr)*((1.0 - eru2)*v_B_r + eru2*dr*htb.v_B*ds_R_B);
  spy_dx = (eru*emu/dr)*(-v_B_r + dr*htb.v_B*ds_R_B);
  spz_dx = (eru*enu/dr)*(-v_B_r + dr*htb.v_B*ds_R_B);

  pxs_dx  = -spx_dx;
  pxpx_dx = (eru/dr)*(2.0*(1.0 - eru2)*(v_C_r - v_D_r)
	  + eru2*dr*htb.v_C*ds_R_C + (1.0 - eru2)*dr*htb.v_D*ds_R_D);
  pxpy_dx = (emu/dr)*((1.0 - 2.0*eru2)*(v_C_r - v_D_r)
          + eru2*dr*(htb.v_C*ds_R_C - htb.v_D*ds_R_D));
  pxpz_dx = (enu/dr)*((1.0 - 2.0*eru2)*(v_C_r - v_D_r)
	  + eru2*dr*(htb.v_C*ds_R_C - htb.v_D*ds_R_D));

  pys_dx  = -spy_dx;
  pypx_dx = pxpy_dx;
  pypy_dx = (eru/dr)*(-2.0*emu2*(v_C_r - v_D_r)
	  + emu2*dr*htb.v_C*ds_R_C + (1.0 - emu2)*dr*htb.v_D*ds_R_D);
  pypz_dx = (eru*emu*enu/dr)*(-2.0*(v_C_r - v_D_r)
          + dr*(htb.v_C*ds_R_C - htb.v_D*ds_R_D));

  pzs_dx  = -spz_dx;
  pzpx_dx = pxpz_dx;
  pzpy_dx = pypz_dx;
  pzpz_dx = (eru/dr)*(-2.0*enu2*(v_C_r - v_D_r)
          + enu2*dr*htb.v_C*ds_R_C + (1.0 - enu2)*dr*htb.v_D*ds_R_D);

/* Fy = [d(T_ab)/dy] */
  ss_dy  = emu*htb.v_A*ds_R_A;
  spx_dy = spy_dx;
  spy_dy = (1.0/dr)*((1.0 - emu2)*v_B_r + emu2*dr*htb.v_B*ds_R_B);
  spz_dy = (emu*enu/dr)*(-v_B_r + dr*htb.v_B*ds_R_B);

  pxs_dy  = -spx_dy;
  pxpx_dy = (emu/dr)*(-2.0*eru2*(v_C_r - v_D_r)
          + eru2*dr*htb.v_C*ds_R_C + (1.0 - eru2)*dr*htb.v_D*ds_R_D);
  pxpy_dy = (eru/dr)*((1.0 - 2.0*emu2)*(v_C_r - v_D_r)
          + emu2*dr*(htb.v_C*ds_R_C - htb.v_D*ds_R_D));
  pxpz_dy = pypz_dx;

  pys_dy  = -spy_dy;
  pypx_dy = pxpy_dy;
  pypy_dy = (emu/dr)*(2.0*(1.0 - emu2)*(v_C_r - v_D_r)
          + emu2*dr*htb.v_C*ds_R_C + (1.0 - emu2)*dr*htb.v_D*ds_R_D);
  pypz_dy = (enu/dr)*((1.0 - 2.0*emu2)*(v_C_r - v_D_r)
          + emu2*dr*(htb.v_C*ds_R_C - htb.v_D*ds_R_D));

  pzs_dy  = -spz_dy;
  pzpx_dy = pxpz_dy;
  pzpy_dy = pypz_dy;
  pzpz_dy = (emu/dr)*(-2.0*enu2*(v_C_r - v_D_r)
          + enu2*dr*htb.v_C*ds_R_C + (1.0 - enu2)*dr*htb.v_D*ds_R_D);

/* Fz = [d(T_ab)/dz] */
  ss_dz  = enu*htb.v_A*ds_R_A;
  spx_dz = (eru*enu/dr)*(-v_B_r + dr*htb.v_B*ds_R_B);
  spy_dz = (emu*enu/dr)*(-v_B_r + dr*htb.v_B*ds_R_B);
  spz_dz = (1.0/dr)*((1.0-enu2)*v_B_r + enu2*dr*htb.v_B*ds_R_B);

  pxs_dz  = -spx_dz;
  pxpx_dz = (enu/dr)*(-2.0*eru2*(v_C_r - v_D_r)
          + eru2*dr*htb.v_C*ds_R_C + (1.0 - eru2)*dr*htb.v_D*ds_R_D);
  pxpy_dz = pypz_dx;
  pxpz_dz = (eru/dr)*((1.0-2.0*enu2)*(v_C_r - v_D_r)
          + enu2*dr*(htb.v_C*ds_R_C - htb.v_D*ds_R_D));

  pys_dz  = -spy_dz;
  pypx_dz = pxpy_dz;
  pypy_dz = (enu/dr)*(-2.0*emu2*(v_C_r - v_D_r)
          + emu2*dr*htb.v_C*ds_R_C + (1.0 - emu2)*dr*htb.v_D*ds_R_D);
  pypz_dz = (emu/dr)*((1.0 - 2.0*enu2)*(v_C_r - v_D_r)
          + enu2*dr*(htb.v_C*ds_R_C - htb.v_D*ds_R_D));

  pzs_dz  = -spz_dz;
  pzpx_dz = pxpz_dz;
  pzpy_dz = pypz_dz;
  pzpz_dz = (enu/dr)*(2.0*(1.0 - enu2)*(v_C_r - v_D_r)
          + enu2*dr*htb.v_C*ds_R_C + (1.0 - enu2)*dr*htb.v_D*ds_R_D);

/* Fx = d(Htb)/dx */
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*0+0] = ss_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*0+1] = spx_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*0+2] = spy_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*0+3] = spz_dx;

  htb.mat2x[16*sys.N*i+4*j+4*sys.N*1+0] = pxs_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*1+1] = pxpx_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*1+2] = pxpy_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*1+3] = pxpz_dx;

  htb.mat2x[16*sys.N*i+4*j+4*sys.N*2+0] = pys_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*2+1] = pypx_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*2+2] = pypy_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*2+3] = pypz_dx;

  htb.mat2x[16*sys.N*i+4*j+4*sys.N*3+0] = pzs_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*3+1] = pzpx_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*3+2] = pzpy_dx;
  htb.mat2x[16*sys.N*i+4*j+4*sys.N*3+3] = pzpz_dx;

/* Fy = d(Htb)/dy */
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*0+0] = ss_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*0+1] = spx_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*0+2] = spy_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*0+3] = spz_dy;

  htb.mat2y[16*sys.N*i+4*j+4*sys.N*1+0] = pxs_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*1+1] = pxpx_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*1+2] = pxpy_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*1+3] = pxpz_dy;

  htb.mat2y[16*sys.N*i+4*j+4*sys.N*2+0] = pys_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*2+1] = pypx_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*2+2] = pypy_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*2+3] = pypz_dy;

  htb.mat2y[16*sys.N*i+4*j+4*sys.N*3+0] = pzs_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*3+1] = pzpx_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*3+2] = pzpy_dy;
  htb.mat2y[16*sys.N*i+4*j+4*sys.N*3+3] = pzpz_dy;

/* Fz = d(Htb)/dz */
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*0+0] = ss_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*0+1] = spx_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*0+2] = spy_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*0+3] = spz_dz;

  htb.mat2z[16*sys.N*i+4*j+4*sys.N*1+0] = pxs_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*1+1] = pxpx_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*1+2] = pxpy_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*1+3] = pxpz_dz;

  htb.mat2z[16*sys.N*i+4*j+4*sys.N*2+0] = pys_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*2+1] = pypx_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*2+2] = pypy_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*2+3] = pypz_dz;

  htb.mat2z[16*sys.N*i+4*j+4*sys.N*3+0] = pzs_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*3+1] = pzpx_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*3+2] = pzpy_dz;
  htb.mat2z[16*sys.N*i+4*j+4*sys.N*3+3] = pzpz_dz;
}

/*--------------------------------------------------------------------------*/
void   diagonal(int i) {

  int ii, jj;

  htb.mat1[16*sys.N*i+4*i+4*sys.N*0+0] = htb.e_s;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*0+1] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*0+2] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*0+3] = 0.0;

  htb.mat1[16*sys.N*i+4*i+4*sys.N*1+0] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*1+1] = htb.e_p;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*1+2] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*1+3] = 0.0;

  htb.mat1[16*sys.N*i+4*i+4*sys.N*2+0] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*2+1] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*2+2] = htb.e_p;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*2+3] = 0.0;

  htb.mat1[16*sys.N*i+4*i+4*sys.N*3+0] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*3+1] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*3+2] = 0.0;
  htb.mat1[16*sys.N*i+4*i+4*sys.N*3+3] = htb.e_p;

  for(ii=0;ii<4;ii++) {
    for(jj=0;jj<4;jj++) {
      htb.mat2x[16*sys.N*i+4*i+4*sys.N*ii+jj] = 0.0;
      htb.mat2y[16*sys.N*i+4*i+4*sys.N*ii+jj] = 0.0;
      htb.mat2z[16*sys.N*i+4*i+4*sys.N*ii+jj] = 0.0;
    }
  }
}

