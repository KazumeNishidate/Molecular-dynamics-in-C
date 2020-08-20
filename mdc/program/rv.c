#include  <stdlib.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/***
*  Calculation of force:
*   Coulomb term calculation using EWALD method will be automatically 
*   SKIPPED when "alpha parameter" sys.a1 is set to 0.0. See also 
*   "manual.IEMD", "control.c" and "rcprcl.c".
***/
void    calc_foc(void)
{
  clear_foc(); 
  real_space();
  reciprocal_space();  
  reciprocal_space3();
}

/*** 
*  Verlet's velocity form: 
*  - L. Verlet, Phys. Rev. 136. A405 (1964).
*  - "Computer Simulation Methods in Theoretical Physics"
*     D. W. Heermann, Springer-Verlag (1989).
***/
void    next_rv_verlet(void)
{
  calc_foc();  /* calculation of force */
  next_v();
  next_r();
}

void	next_v(void) /* Verlet's velocity form */ 
{
  int i;
  double m2;

  for(i=0; i<sys.N; i++) {
    m2 = 2.0*ion_m(i); 
    sys.vx[i] += sys.dt*(sys.fx[i]+sys.fx0[i])/m2;
    sys.vy[i] += sys.dt*(sys.fy[i]+sys.fy0[i])/m2; 
    sys.vz[i] += sys.dt*(sys.fz[i]+sys.fz0[i])/m2; 
  }
}

void	next_r(void) /* Verlet's velocity form */ 
{
  int i;
  double m2;

  for(i=0; i<sys.N; i++) {

    m2 = 2.0*ion_m(i); 
    sys.rx[i] += sys.dt * sys.vx[i] + sys.dt2*sys.fx[i]/m2;
    sys.ry[i] += sys.dt * sys.vy[i] + sys.dt2*sys.fy[i]/m2;
    sys.rz[i] += sys.dt * sys.vz[i] + sys.dt2*sys.fz[i]/m2;

    /* cyclic boundary condition */
    if(sys.rx[i] > sys.Lx) sys.rx[i] = sys.rx[i]-sys.Lx;
    if(sys.ry[i] > sys.Ly) sys.ry[i] = sys.ry[i]-sys.Ly;
    if(sys.rz[i] > sys.Lz) sys.rz[i] = sys.rz[i]-sys.Lz;
    if(sys.rx[i] < 0.0) sys.rx[i] = sys.Lx+sys.rx[i];
    if(sys.ry[i] < 0.0) sys.ry[i] = sys.Ly+sys.ry[i];
    if(sys.rz[i] < 0.0) sys.rz[i] = sys.Lz+sys.rz[i];
  }     
}

void   	calc_kin(void) /* calculation of kinetic energy */
{ 
  int i;
  double kx,ky,kz,m;

  kx=ky=kz=0.0;

  for(i=0; i<sys.N; i++){
    m = ion_m(i);
    kx += m*sys.vx[i]*sys.vx[i];
    ky += m*sys.vy[i]*sys.vy[i];
    kz += m*sys.vz[i]*sys.vz[i];
  }

  sys.kinX = kx /2.0 ;
  sys.kinY = ky /2.0 ;
  sys.kinZ = kz /2.0 ;
  sys.kin = sys.kinX+sys.kinY+sys.kinZ;
}

void    moment_correction(void)
{
  int i;
  double vx_sum, vy_sum, vz_sum;
  double mean_vx, mean_vy, mean_vz;
  double m;

  vx_sum=0.0;
  vy_sum=0.0;
  vz_sum=0.0;

  for(i=0; i<sys.N; i++){
    m=ion_m(i);
    vx_sum += m*sys.vx[i];
    vy_sum += m*sys.vy[i];
    vz_sum += m*sys.vz[i];
  }

  mean_vx = vx_sum/((double)sys.N);
  mean_vy = vy_sum/((double)sys.N);
  mean_vz = vz_sum/((double)sys.N);

  for(i=0; i<sys.N; i++){
    m=ion_m(i);
    sys.vx[i] -= mean_vx/m;
    sys.vy[i] -= mean_vy/m;
    sys.vz[i] -= mean_vz/m;
  }
}

/***
*  Gear's 7th-order F representation:
*  - G. Ciccotti and W. G. Hoover eds, "Molecular Dynamics Simulation
*    of Statistical-Mechanical Systems", North-Holland Physics (1986).
*  - Y. Hiwatari, Solid State Physics (KOTAIBUTSURI), Vol. 24, No. 277, 
*    pp108(242)-118(252) (1989). 
***/
void next_rv_gear(void) 
{
  static int cnt=0;
  static double A13, A14, A15, A16, A17;
  static double A23, A24, A25, A26, A27;
  static double A33, A34, A35, A36, A37;
  static double a1, a2;
  static double dt, dt2, dt22;

  int i;
  double mass;
  double pre_ax, pre_ay, pre_az;

  if(cnt==0) {  /* make a preparation at first time step */
    A13 =  1427.0/720.0;  A14 =  -133.0/60.0;   A15 =   241.0/120.0;  
    A16 =  -173.0/180.0;  A17 =     3.0/16.0;

    A23 =  1901.0/360.0;  A24 = -1387.0/180.0;  A25 =   109.0/15.0;   
    A26 =  -637.0/180.0;  A27 =   251.0/360.0;

    A33 = 5.0;  A34 = -10.0;  A35 = 10.0;  A36 = -5.0;  A37 = 1.0;

    a1 = 863.0/6048.0;   a2 = 665.0/1008.0;
    gear.anx  = (double *)calloc(sys.N, sizeof(double)); 
    gear.any  = (double *)calloc(sys.N, sizeof(double)); 
    gear.anz  = (double *)calloc(sys.N, sizeof(double)); 
    gear.an1x = (double *)calloc(sys.N, sizeof(double)); 
    gear.an1y = (double *)calloc(sys.N, sizeof(double)); 
    gear.an1z = (double *)calloc(sys.N, sizeof(double)); 
    gear.an2x = (double *)calloc(sys.N, sizeof(double)); 
    gear.an2y = (double *)calloc(sys.N, sizeof(double)); 
    gear.an2z = (double *)calloc(sys.N, sizeof(double)); 
    gear.an3x = (double *)calloc(sys.N, sizeof(double)); 
    gear.an3y = (double *)calloc(sys.N, sizeof(double)); 
    gear.an3z = (double *)calloc(sys.N, sizeof(double)); 
    gear.an4x = (double *)calloc(sys.N, sizeof(double)); 
    gear.an4y = (double *)calloc(sys.N, sizeof(double)); 
    gear.an4z = (double *)calloc(sys.N, sizeof(double)); 

    dt = sys.dt;
    dt2  = dt/2.0;
    dt22 = dt*dt/2.0;
    cnt = 1;
  }

  for(i=0; i<sys.N; i++) {   /* predictor */
    sys.rx[i] += dt*sys.vx[i]+
      dt22*( A13*gear.anx[i]  + A14*gear.an1x[i] + A15*gear.an2x[i] +
                                A16*gear.an3x[i] + A17*gear.an4x[i]);
    sys.ry[i] += dt*sys.vy[i]+
      dt22*( A13*gear.any[i]  + A14*gear.an1y[i] + A15*gear.an2y[i] +
                                A16*gear.an3y[i] + A17*gear.an4y[i]);
    sys.rz[i] += dt*sys.vz[i]+
      dt22*( A13*gear.anz[i]  + A14*gear.an1z[i] + A15*gear.an2z[i] +
                                A16*gear.an3z[i] + A17*gear.an4z[i]);

    /* cyclic boundary condition */
    if(sys.rx[i] > sys.Lx) sys.rx[i] = sys.rx[i]-sys.Lx;
    if(sys.ry[i] > sys.Ly) sys.ry[i] = sys.ry[i]-sys.Ly;
    if(sys.rz[i] > sys.Lz) sys.rz[i] = sys.rz[i]-sys.Lz;
    if(sys.rx[i] < 0.0) sys.rx[i] = sys.Lx+sys.rx[i];
    if(sys.ry[i] < 0.0) sys.ry[i] = sys.Ly+sys.ry[i];
    if(sys.rz[i] < 0.0) sys.rz[i] = sys.Lz+sys.rz[i];

    sys.vx[i] += dt2*(A23*gear.anx[i] + A24*gear.an1x[i] +
		     A25*gear.an2x[i] + A26*gear.an3x[i] + A27*gear.an4x[i]);
    sys.vy[i] += dt2*(A23*gear.any[i] + A24*gear.an1y[i] +
		     A25*gear.an2y[i] + A26*gear.an3y[i] + A27*gear.an4y[i]);
    sys.vz[i] += dt2*(A23*gear.anz[i] + A24*gear.an1z[i] +
		     A25*gear.an2z[i] + A26*gear.an3z[i] + A27*gear.an4z[i]);
  }

  calc_foc();  /* calculation of force */

  for(i=0; i<sys.N; i++) {   /* corrector */
    mass = ion_m(i);

    pre_ax = A33*gear.anx[i] + A34*gear.an1x[i] + 
      A35*gear.an2x[i] + A36*gear.an3x[i] + A37*gear.an4x[i];

    pre_ay = A33*gear.any[i] + A34*gear.an1y[i] + 
      A35*gear.an2y[i] + A36*gear.an3y[i] + A37*gear.an4y[i];

    pre_az = A33*gear.anz[i] + A34*gear.an1z[i] + 
      A35*gear.an2z[i] + A36*gear.an3z[i] + A37*gear.an4z[i];

    gear.an4x[i] = gear.an3x[i];
    gear.an4y[i] = gear.an3y[i];
    gear.an4z[i] = gear.an3z[i];

    gear.an3x[i] = gear.an2x[i];
    gear.an3y[i] = gear.an2y[i];
    gear.an3z[i] = gear.an2z[i];

    gear.an2x[i] = gear.an1x[i];
    gear.an2y[i] = gear.an1y[i];
    gear.an2z[i] = gear.an1z[i];

    gear.an1x[i] = gear.anx[i];
    gear.an1y[i] = gear.any[i];
    gear.an1z[i] = gear.anz[i];

    gear.anx[i]  = sys.fx[i]/mass;
    gear.any[i]  = sys.fy[i]/mass;
    gear.anz[i]  = sys.fz[i]/mass;

    sys.rx[i] += a1*dt22*(gear.anx[i] - pre_ax);
    sys.ry[i] += a1*dt22*(gear.any[i] - pre_ay);
    sys.rz[i] += a1*dt22*(gear.anz[i] - pre_az);

    /* cyclic boundary condition */
    if(sys.rx[i] > sys.Lx) sys.rx[i] = sys.rx[i]-sys.Lx;
    if(sys.ry[i] > sys.Ly) sys.ry[i] = sys.ry[i]-sys.Ly;
    if(sys.rz[i] > sys.Lz) sys.rz[i] = sys.rz[i]-sys.Lz;
    if(sys.rx[i] < 0.0) sys.rx[i] = sys.Lx+sys.rx[i];
    if(sys.ry[i] < 0.0) sys.ry[i] = sys.Ly+sys.ry[i];
    if(sys.rz[i] < 0.0) sys.rz[i] = sys.Lz+sys.rz[i];

    sys.vx[i] += a2*dt2*(gear.anx[i] - pre_ax);
    sys.vy[i] += a2*dt2*(gear.any[i] - pre_ay);
    sys.vz[i] += a2*dt2*(gear.anz[i] - pre_az);
  }
}



