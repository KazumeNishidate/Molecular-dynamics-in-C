#include  <stdlib.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

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

/***
*  Calculation of force:
*   Coulomb term calculation using EWALD method will be automatically 
*   SKIPPED when "ctl.skip_ewald" is set to 1 in "control.c".
*        calc Ewald : ctl.skip_ewald = 0
*        skip Ewald : ctl.skip_ewald = 1
***/
void    calc_foc(void)
{
  clear_foc(); 

  calc_sigma();  /* sigma matrix generation */
  set_cutoff();  /* cutoff radius  */ 

  real_space();  /* repulsive part */

  cell_force();

  if(ctl.skip_ewald == 1) return;  /* 0:calc Ewald, 1:skip Ewald */

  ewald1();      /*                1st term */
  ewald2();      /*                2nd term */
  ewald3();      /*                3rd term */
}

void	next_v(void) /* Verlet's velocity form */ 
{
  int i;
  double m2;

  for(i=0; i<sys.N; i++) {
    m2 = 2.0*ion_m(i); 
    cell.vx[i] += sys.dt*(cell.fx[i]+cell.fx0[i])/m2;
    cell.vy[i] += sys.dt*(cell.fy[i]+cell.fy0[i])/m2; 
    cell.vz[i] += sys.dt*(cell.fz[i]+cell.fz0[i])/m2; 
  }
}

void	next_r(void) /* Verlet's velocity form */ 
{
  int i;
  double m2;

  for(i=0; i<sys.N; i++) {
    m2 = 2.0*ion_m(i); 
    cell.rx[i] += sys.dt * cell.vx[i] + sys.dt2*cell.fx[i]/m2;
    cell.ry[i] += sys.dt * cell.vy[i] + sys.dt2*cell.fy[i]/m2;
    cell.rz[i] += sys.dt * cell.vz[i] + sys.dt2*cell.fz[i]/m2;
    cyclic_bc(i);    /* cyclic boundary condition */
  }     
}

void   	calc_kin(void) /* calculation of kinetic energy tensor [mv^2] */
{ 
  int i;
  double kxx, kxy, kxz;
  double kyx, kyy, kyz;
  double kzx, kzy, kzz;
  double m;

  kxx = 0.0;   kxy = 0.0;   kxz = 0.0; 
  kyx = 0.0;   kyy = 0.0;   kyz = 0.0; 
  kzx = 0.0;   kzy = 0.0;   kzz = 0.0; 

  for(i=0; i<sys.N; i++){
    m = ion_m(i);
    kxx += m*cell.vx[i]*cell.vx[i];
    kxy += m*cell.vx[i]*cell.vy[i];
    kxz += m*cell.vx[i]*cell.vz[i];

    kyy += m*cell.vy[i]*cell.vy[i];
    kyz += m*cell.vy[i]*cell.vz[i];

    kzz += m*cell.vz[i]*cell.vz[i];
  }

  kyx = kxy;   kzx = kxz;  kzy = kyz;

  pt.kin_tensor[0][0] = kxx;
  pt.kin_tensor[0][1] = kxy;
  pt.kin_tensor[0][2] = kxz;

  pt.kin_tensor[1][0] = kyx;
  pt.kin_tensor[1][1] = kyy;
  pt.kin_tensor[1][2] = kyz;

  pt.kin_tensor[2][0] = kzx;
  pt.kin_tensor[2][1] = kzy;
  pt.kin_tensor[2][2] = kzz;

  sys.kin = (kxx + kyy + kzz)/2.0;
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
    vx_sum += m*cell.vx[i];
    vy_sum += m*cell.vy[i];
    vz_sum += m*cell.vz[i];
  }

  mean_vx = vx_sum/((double)sys.N);
  mean_vy = vy_sum/((double)sys.N);
  mean_vz = vz_sum/((double)sys.N);

  for(i=0; i<sys.N; i++){
    m=ion_m(i);
    cell.vx[i] -= mean_vx/m;
    cell.vy[i] -= mean_vy/m;
    cell.vz[i] -= mean_vz/m;
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
    cell.rx[i] += dt*cell.vx[i]+
      dt22*( A13*gear.anx[i]  + A14*gear.an1x[i] + A15*gear.an2x[i] +
                                A16*gear.an3x[i] + A17*gear.an4x[i]);
    cell.ry[i] += dt*cell.vy[i]+
      dt22*( A13*gear.any[i]  + A14*gear.an1y[i] + A15*gear.an2y[i] +
                                A16*gear.an3y[i] + A17*gear.an4y[i]);
    cell.rz[i] += dt*cell.vz[i]+
      dt22*( A13*gear.anz[i]  + A14*gear.an1z[i] + A15*gear.an2z[i] +
                                A16*gear.an3z[i] + A17*gear.an4z[i]);

    cyclic_bc(i); /* cyclic boundary condition */

    cell.vx[i] += dt2*(A23*gear.anx[i] + A24*gear.an1x[i] +
		     A25*gear.an2x[i] + A26*gear.an3x[i] + A27*gear.an4x[i]);
    cell.vy[i] += dt2*(A23*gear.any[i] + A24*gear.an1y[i] +
		     A25*gear.an2y[i] + A26*gear.an3y[i] + A27*gear.an4y[i]);
    cell.vz[i] += dt2*(A23*gear.anz[i] + A24*gear.an1z[i] +
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

    gear.anx[i]  = cell.fx[i]/mass;
    gear.any[i]  = cell.fy[i]/mass;
    gear.anz[i]  = cell.fz[i]/mass;

    cell.rx[i] += a1*dt22*(gear.anx[i] - pre_ax);
    cell.ry[i] += a1*dt22*(gear.any[i] - pre_ay);
    cell.rz[i] += a1*dt22*(gear.anz[i] - pre_az);

    cyclic_bc(i); /* cyclic boundary condition */

    cell.vx[i] += a2*dt2*(gear.anx[i] - pre_ax);
    cell.vy[i] += a2*dt2*(gear.any[i] - pre_ay);
    cell.vz[i] += a2*dt2*(gear.anz[i] - pre_az);
  }
}

