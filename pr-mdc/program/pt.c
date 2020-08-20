#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*--------------------------------------------------------------------
*  Calculation of pressure tensor P_xx, P_yy, and P_zz using VIRIAL
*---------------------------------------------------------------------*/
void   calc_press(void)
{
  double **pp;
  int i, j;

  pp = pt.pres_tensor; 
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      pp[i][j] = (2.0*pt.kin_tensor[i][j] + pt.virial[i][j])/cell.vol;
      pp[i][j] = (2.0*pt.kin_tensor[i][j] + pt.virial[i][j])/cell.vol;
      pp[i][j] = (2.0*pt.kin_tensor[i][j] + pt.virial[i][j])/cell.vol;
    }
  }

  /* current pressure [ = Trace[Pab]/3 ] */
  pt.pres_trace =  ( pp[0][0] + pp[1][1] + pp[2][2])/3.0;

}

void   control_press_PR() {
  int i, j, k;
  double **sg, **f, **f0;
  double ff[3][3];

  sg = cell.sigma;
  f  = pr.f;
  f0 = pr.f0;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++) /* pi-P */
      ff[i][j] = pt.pres_tensor[i][j] - pt.pres_set[i][j];

  for(i=0;i<3;i++)
    for(j=0;j<3;j++) {
      pr.f0[i][j] = pr.f[i][j];
      pr.f[i][j]  = 0.0;
      for(k=0;k<3;k++) /* F = Wh'' = (pi-P).SIGMA */  
	pr.f[i][j] += ff[i][k]*sg[k][j];
    }

  for(i=0;i<3;i++) 
    for(j=0;j<3;j++)
      pr.v[i][j] += pr.dt*(pr.f[i][j]+pr.f0[i][j])/(2.0*pr.mass);

  for(i=0;i<3;i++) 
    for(j=0;j<3;j++)
      cell.hm[i][j] += 
	pr.dt * pr.v[i][j] + (pr.dt*pr.dt)*pr.f[i][j]/(2.0*pr.mass);

  calc_cell(); /* evaluate the cell matrix : Hm -> Hmr */

  for(i=0; i<sys.N; i++) {
    cyclic_bc(i);    /* cyclic boundary condition */
  }     

}

/*****
*   pressure control by forced scaling method
*    K. Kawamura, "Molecular Dynamics Simulations" F. Yonezawa ed.,
*    Solid-State Sciences Vol. 103, Springer-Verlag, p.88, (1992).
*****/
void   control_press(int d_step)
{
  static double oldPx = 0.0, oldPy = 0.0, oldPz = 0.0;
  double averaged_Px, averaged_Py, averaged_Pz;
  double target_Px, target_Py, target_Pz;
  static short count = 0;
  double SX, SY, SZ;  /* scale factor */
  int i;

  count++;

  oldPx += pt.pres_tensor[0][0];
  oldPy += pt.pres_tensor[1][1];
  oldPz += pt.pres_tensor[2][2];
	
  if(sys.step % d_step == 0) {

    averaged_Px = oldPx / ((double)count);
    averaged_Py = oldPy / ((double)count);
    averaged_Pz = oldPz / ((double)count);
    oldPx = 0.0;
    oldPy = 0.0;
    oldPz = 0.0;
    count = 0;

    target_Px = pt.pres_set[0][0];
    target_Py = pt.pres_set[1][1];
    target_Pz = pt.pres_set[2][2];

    /* evaluate a cell size scaling factor [SX, SY, SZ]        */
    /* where "pt.pres_scalling_factor" is a constant (=0.5).   */
    SX = 1.0 + atan(averaged_Px - target_Px)*pt.pres_scalling_factor;
    SY = 1.0 + atan(averaged_Py - target_Py)*pt.pres_scalling_factor;
    SZ = 1.0 + atan(averaged_Pz - target_Pz)*pt.pres_scalling_factor;

    /* scaling */
    cell.hm[0][0] *= SX;
    cell.hm[1][1] *= SY;
    cell.hm[2][2] *= SZ;

    /* position shift after the cell size scaling */
    for(i=0; i<sys.N; i++){
      cell.rx[i] *= SX;
      cell.ry[i] *= SY;
      cell.rz[i] *= SZ;
    }

   calc_cell(); /* evaluate the cell matrix : Hm -> Hmr        */

   printf("---------------------------------------------------------------\n");
   printf(">> Pressure scaling\n");
  }
}

void   control_temp(int d_step, double temp)
{
  double scale, tempK, averaged_K = 0.0;
  static double oldk = 0.0;    
  static short count = 0;
  int i;

  count++;
  oldk += sys.kin;

  if(sys.step % d_step == 0) { /* evaluate the averaged kinetic energy */
    averaged_K = oldk / (double) count; 
    oldk = 0.0;
    count = 0;

    tempK  = temp/sys.e2t;                /* temperature -> [energy]     */
    scale = sqrt(tempK / averaged_K);     /* evaluate the scaling factor */
    for(i=0; i< sys.N; i++) {
      cell.vx[i] *= scale;
      cell.vy[i] *= scale;
      cell.vz[i] *= scale;
    }

    /*   printf("---------------------------------------------------------------\n"); */
    /*   printf(">> Temperature scaling\n"); */

  }
}

