#include  <stdlib.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*****
*    Parrinello Rahman method
*    M. Parrinello and A. Rahman, Phys. Rev. Lett., Vo.45, 1196 (1980).
*    M. Parrinello and A. Rahman, J. Appl. Phys., Vo.52, 7182 (1981).
*
*    # integrate the cell EOM using verlet velocity form. 
******/
void evaluate_PR_param(void) {
  int i;
  double mass = 0.0;

  for(i=0;i<sys.N;i++) {
    mass += ion_m(i);
  }
  pr.mass = 3.0*mass/(4.0*PIPI); /* cell mass */
  pr.dt = sys.dt;                /* cell dt   */
}

/*----------------------------------------------------------------*/
/* cutoff length = the radius of the sphere in MD super cell box  */
/*----------------------------------------------------------------*/
void set_cutoff(void) {
  double h1, h2, h3;
  double min=0.0;
  double **g;

  g = cell.g;

  /* length between the opposite MD cell faces */
  h1 = 1.0/sqrt(g[0][0]*g[0][0]+g[0][1]*g[0][1]+g[0][2]*g[0][2]);
  h2 = 1.0/sqrt(g[1][0]*g[1][0]+g[1][1]*g[1][1]+g[1][2]*g[1][2]);
  h3 = 1.0/sqrt(g[2][0]*g[2][0]+g[2][1]*g[2][1]+g[2][2]*g[2][2]);

  if((h1 <= h2) && (h1 <= h3)) min = h1;
  if((h2 <= h1) && (h2 <= h3)) min = h2;
  if((h3 <= h1) && (h3 <= h2)) min = h3;

  if(sys.step>10) return;

  sys.cutoff  = min/2.0;  /* real space cutoff [A] */
  sys.cutoff2 = sys.cutoff*sys.cutoff;

  if(ctl.skip_ewald == 1) return; /* 0:calc Ewald, 1:skip Ewald */

  sys.a1  = 3.8/min;  /* note that "alpha!=5.6/min" in general */
}

/*-----------------------------------------------*/
/* cell force calculation                        */
/* (d/dt) p = F - [h^T]^(-1) * [(d/dt)h]^T * p   */
/*-----------------------------------------------*/
void cell_force(void) {
  int i, j, k;
  double **Hm, **HmT, **HmTr;
  double **Hmv, **HmvT, **HmF;
  double det;
  double mass;

  Hm      = cell.hm;
  HmT     = cell.hmt;
  HmTr    = cell.hmtr;
  Hmv     = pr.v;
  HmvT    = cell.hmvt; 
  HmF     = cell.hmf;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++) /* transvers */
      HmT[i][j] = Hm[j][i];

  det = HmT[0][0]*(HmT[1][1]*HmT[2][2] - HmT[2][1]*HmT[1][2])
    + HmT[1][0]*(HmT[2][1]*HmT[0][2] - HmT[0][1]*HmT[2][2])
    + HmT[2][0]*(HmT[0][1]*HmT[1][2] - HmT[1][1]*HmT[0][2]);

  HmTr[0][0] = (-(HmT[1][2]*HmT[2][1]) + HmT[1][1]*HmT[2][2])/det;
  HmTr[1][0] = (  HmT[1][2]*HmT[2][0]  - HmT[1][0]*HmT[2][2])/det;
  HmTr[2][0] = (-(HmT[1][1]*HmT[2][0]) + HmT[1][0]*HmT[2][1])/det;
  HmTr[0][1] = (  HmT[0][2]*HmT[2][1]  - HmT[0][1]*HmT[2][2])/det;
  HmTr[1][1] = (-(HmT[0][2]*HmT[2][0]) + HmT[0][0]*HmT[2][2])/det;
  HmTr[2][1] = (  HmT[0][1]*HmT[2][0]  - HmT[0][0]*HmT[2][1])/det;
  HmTr[0][2] = (-(HmT[0][2]*HmT[1][1]) + HmT[0][1]*HmT[1][2])/det;
  HmTr[1][2] = (  HmT[0][2]*HmT[1][0]  - HmT[0][0]*HmT[1][2])/det;
  HmTr[2][2] = (-(HmT[0][1]*HmT[1][0]) + HmT[0][0]*HmT[1][1])/det;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      HmvT[i][j] = Hmv[j][i];

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      HmF[i][j] = 0.0;
      
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for(k=0;k<3;k++) /* [h^T]^(-1) * [(d/dt)h]^T */
	HmF[i][j] += HmTr[i][k]*HmvT[k][j];

  for(i=0;i<sys.N;i++) {
    mass = ion_m(i);

    cell.fx[i] -= mass*(HmF[0][0]*cell.vx[i] + HmF[1][0]*cell.vy[i]
                                       + HmF[2][0]*cell.vz[i]);

    cell.fy[i] -= mass*(HmF[0][1]*cell.vx[i] + HmF[1][1]*cell.vy[i]
                                       + HmF[2][1]*cell.vz[i]);

    cell.fz[i] -= mass*(HmF[0][2]*cell.vx[i] + HmF[1][2]*cell.vy[i]
                                       + HmF[2][2]*cell.vz[i]);
  }
}

/*---------------------------------------------------*/
/* calculation of the sigma matrix                   */
/*                                                   */
/* Cartesian coordinate{r} <-> scaled coordinate {s} */
/*                                                   */
/* H = {a, b, c}, r = H*s                            */
/*   r : Cartesian coordinates                       */
/*   s : reduced coordinates                         */
/*   H : MD super cell matrix                        */
/*                                                   */
/*  {r} = H * {s}                                    */
/*  {s} = H^(-1) * {r}                               */
/*                                                   */
/*  |rx|   |h00 h10 h20| |sx|                        */
/*  |ry| = |h01 h11 h21|*|sy|                        */
/*  |rz|   |h02 h12 h22| |sz|                        */
/*                                                   */
/*                                                   */
/* sigma ={ b x c, c x a, a x b }                    */
/* g     = sigma/(cell volume)                       */
/*                                                   */
/* sigma describes the size and orientation of       */
/* the MD super cell faces                           */ 
/*                                                   */
/* reciprocal lattice vector g = {gx,gy,gz}          */
/*                                                   */
/*  |gx|     2Pi    |sig00 sig10 sig20| |h|          */
/*  |gy| = -------- |sig01 sig11 sig21|*|k|          */
/*  |gz|    volume  |sig02 sig12 sig22| |l|          */
/*                                                   */
/*  h,k,l => integer                                 */
/*---------------------------------------------------*/
void calc_sigma(void) { 
  double **Hm, **sg, **g;
  int i, j;

  Hm = cell.hm;
  g  = cell.g;
  sg = cell.sigma;

  /* sigma1 = {b x c} */
  sg[0][0] = (Hm[1][1]*Hm[2][2] - Hm[1][2]*Hm[2][1]);
  sg[0][1] = (Hm[1][2]*Hm[2][0] - Hm[1][0]*Hm[2][2]);
  sg[0][2] = (Hm[1][0]*Hm[2][1] - Hm[1][1]*Hm[2][0]);

  /* sigma2 = {c x a} */
  sg[1][0] = (Hm[2][1]*Hm[0][2] - Hm[2][2]*Hm[0][1]);
  sg[1][1] = (Hm[2][2]*Hm[0][0] - Hm[2][0]*Hm[0][2]);
  sg[1][2] = (Hm[2][0]*Hm[0][1] - Hm[2][1]*Hm[0][0]);

  /* sigma3 = {a x b} */
  sg[2][0] = (Hm[0][1]*Hm[1][2] - Hm[0][2]*Hm[1][1]);
  sg[2][1] = (Hm[0][2]*Hm[1][0] - Hm[0][0]*Hm[1][2]);
  sg[2][2] = (Hm[0][0]*Hm[1][1] - Hm[0][1]*Hm[1][0]);

  for(i=0;i<3;i++) 
    for(j=0;j<3;j++) /* {g} = sigma/(cell volume) */
      g[i][j] = sg[i][j]/cell.vol;

}

/*---------------------------------------------------*/
/* Cartesian coordinate{r} <-> scaled coordinate {s} */
/*                                                   */
/*  {r} = H * {s}                                    */
/*  {s} = H^(-1) * {r}                               */
/*                                                   */
/*  |rx|   |h00 h10 h20| |sx|                        */
/*  |ry| = |h01 h11 h21|*|sy|                        */
/*  |rz|   |h02 h12 h22| |sz|                        */
/*                                                   */
/*  |sx|     1    |rh00 rh10 rh20| |rx|              */
/*  |sy| = ------*|rh01 rh11 rh21|*|ry|              */
/*  |sz|   det(H) |rh02 rh12 rh22| |rz|              */
/*                                                   */
/*    det(H) = h00*(h11*h22 - h21*h12)               */
/*           + h10*(h21*h02 - h01*h22)               */
/*           + h20*(h01*h12 - h11*h02)               */
/*                                                   */
/*      rh00 = -(h12*h21) + h11*h22                  */
/*      rh10 =   h12*h20  - h10*h22                  */
/*      rh20 = -(h11*h20) + h10*h21                  */
/*      rh01 =   h02*h21  - h01*h22                  */
/*      rh11 = -(h02*h20) + h00*h22                  */
/*      rh21 =   h01*h20  - h00*h21                  */
/*      rh02 = -(h02*h11) + h01*h12                  */
/*      rh12 =   h02*h10  - h00*h12                  */
/*      rh22 = -(h01*h10) + h00*h11                  */
/*---------------------------------------------------*/
void calc_cell(void) {
  double det;
  double **Hm, **Hmr;
  double a_b, b_c, c_a;
  double cos_ab, cos_bc, cos_ca;

  Hm  = cell.hm;
  Hmr = cell.hmr;

  det = Hm[0][0]*(Hm[1][1]*Hm[2][2] - Hm[2][1]*Hm[1][2])
    + Hm[1][0]*(Hm[2][1]*Hm[0][2] - Hm[0][1]*Hm[2][2])
    + Hm[2][0]*(Hm[0][1]*Hm[1][2] - Hm[1][1]*Hm[0][2]);

  if(det == 0.0) {
    printf("cell matrix determinant = 0, exiting\n");
    exit(0);
  }

  cell.vol = fabs(det);

  /* evaluate the cell matrix : Hm -> Hmr */
  Hmr[0][0] = (-(Hm[1][2]*Hm[2][1]) + Hm[1][1]*Hm[2][2])/det;
  Hmr[1][0] = (  Hm[1][2]*Hm[2][0]  - Hm[1][0]*Hm[2][2])/det;
  Hmr[2][0] = (-(Hm[1][1]*Hm[2][0]) + Hm[1][0]*Hm[2][1])/det;
  Hmr[0][1] = (  Hm[0][2]*Hm[2][1]  - Hm[0][1]*Hm[2][2])/det;
  Hmr[1][1] = (-(Hm[0][2]*Hm[2][0]) + Hm[0][0]*Hm[2][2])/det;
  Hmr[2][1] = (  Hm[0][1]*Hm[2][0]  - Hm[0][0]*Hm[2][1])/det;
  Hmr[0][2] = (-(Hm[0][2]*Hm[1][1]) + Hm[0][1]*Hm[1][2])/det;
  Hmr[1][2] = (  Hm[0][2]*Hm[1][0]  - Hm[0][0]*Hm[1][2])/det;
  Hmr[2][2] = (-(Hm[0][1]*Hm[1][0]) + Hm[0][0]*Hm[1][1])/det;

  cell.Lx = sqrt(Hm[0][0]*Hm[0][0]+Hm[0][1]*Hm[0][1]+Hm[0][2]*Hm[0][2]);
  cell.Ly = sqrt(Hm[1][0]*Hm[1][0]+Hm[1][1]*Hm[1][1]+Hm[1][2]*Hm[1][2]);
  cell.Lz = sqrt(Hm[2][0]*Hm[2][0]+Hm[2][1]*Hm[2][1]+Hm[2][2]*Hm[2][2]);

  a_b = Hm[0][0]*Hm[1][0]+Hm[0][1]*Hm[1][1]+Hm[0][2]*Hm[1][2];
  b_c = Hm[1][0]*Hm[2][0]+Hm[1][1]*Hm[2][1]+Hm[1][2]*Hm[2][2];
  c_a = Hm[2][0]*Hm[0][0]+Hm[2][1]*Hm[0][1]+Hm[2][2]*Hm[0][2];

  cos_ab = a_b/(cell.Lx*cell.Ly);
  cos_bc = b_c/(cell.Ly*cell.Lz);
  cos_ca = c_a/(cell.Lz*cell.Lx);

  cell.angle_ab = RADIAN_TO_DEGREE(acos(cos_ab));
  cell.angle_bc = RADIAN_TO_DEGREE(acos(cos_bc));
  cell.angle_ca = RADIAN_TO_DEGREE(acos(cos_ca));

}

void cyclic_bc(int i)
{
  double *sx, *sy, *sz;

  sx = cell.sx;  sy = cell.sy;  sz = cell.sz;

  r2s(i); /* Cartesian coordinate -> scaled coordinate */

  if(sx[i] >= 1.0) sx[i] -= 1.0;   /* cyclic boundary condition */
  if(sy[i] >= 1.0) sy[i] -= 1.0;
  if(sz[i] >= 1.0) sz[i] -= 1.0;
  if(sx[i] < 0.0)  sx[i] += 1.0;
  if(sy[i] < 0.0)  sy[i] += 1.0;
  if(sz[i] < 0.0)  sz[i] += 1.0;

  s2r(i); /* scaled coordinate -> Cartesian coordinate */

}

void r2s(int i) { /* Cartesian coordinate -> scaled coordinate */
  double **Hm, **Hmr;
  double *rx, *ry, *rz;
  double *sx, *sy, *sz;

  Hm  = cell.hm;  Hmr = cell.hmr;

  rx = cell.rx;  ry = cell.ry;  rz = cell.rz;
  sx = cell.sx;  sy = cell.sy;  sz = cell.sz;

  sx[i] = Hmr[0][0]*rx[i] + Hmr[1][0]*ry[i] + Hmr[2][0]*rz[i];
  sy[i] = Hmr[0][1]*rx[i] + Hmr[1][1]*ry[i] + Hmr[2][1]*rz[i];
  sz[i] = Hmr[0][2]*rx[i] + Hmr[1][2]*ry[i] + Hmr[2][2]*rz[i];

}

void s2r(int i) { /* scaled coordinate -> Cartesian coordinate */
  double **Hm;
  double *rx, *ry, *rz;
  double *sx, *sy, *sz;

  Hm  = cell.hm;

  rx = cell.rx;  ry = cell.ry;  rz = cell.rz;
  sx = cell.sx;  sy = cell.sy;  sz = cell.sz;

  rx[i] = Hm[0][0]*sx[i] + Hm[1][0]*sy[i] + Hm[2][0]*sz[i];
  ry[i] = Hm[0][1]*sx[i] + Hm[1][1]*sy[i] + Hm[2][1]*sz[i];
  rz[i] = Hm[0][2]*sx[i] + Hm[1][2]*sy[i] + Hm[2][2]*sz[i];

}

void cyclic_bc_dr(double *dx, double *dy, double *dz)
{
  double sdx, sdy, sdz;

  double **Hm, **Hmr;

  Hm  = cell.hm;  Hmr = cell.hmr;

  /* Cartesian coordinate -> scaled coordinate */
  sdx = Hmr[0][0]*(*dx) + Hmr[1][0]*(*dy) + Hmr[2][0]*(*dz);
  sdy = Hmr[0][1]*(*dx) + Hmr[1][1]*(*dy) + Hmr[2][1]*(*dz);
  sdz = Hmr[0][2]*(*dx) + Hmr[1][2]*(*dy) + Hmr[2][2]*(*dz);

  if(sdx >=  0.5) sdx -= 1.0; /* cyclic BC in scaled cell */
  if(sdy >=  0.5) sdy -= 1.0;
  if(sdz >=  0.5) sdz -= 1.0;
  if(sdx <  -0.5) sdx += 1.0;
  if(sdy <  -0.5) sdy += 1.0;
  if(sdz <  -0.5) sdz += 1.0;

  /* scaled coordinate -> Cartesian coordinate */
  *dx = Hm[0][0]*sdx + Hm[1][0]*sdy + Hm[2][0]*sdz;
  *dy = Hm[0][1]*sdx + Hm[1][1]*sdy + Hm[2][1]*sdz;
  *dz = Hm[0][2]*sdx + Hm[1][2]*sdy + Hm[2][2]*sdz;
}

