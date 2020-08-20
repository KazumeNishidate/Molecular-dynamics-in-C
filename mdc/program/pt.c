#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*--------------------------------------------------------------------
*  Calculation of pressure tensor P_xx, P_yy, and P_zz using VIRIAL
*---------------------------------------------------------------------*/
void   calc_press(void)
{
  double md_cell_volume;

  md_cell_volume = sys.Lx*sys.Ly*sys.Lz;

  sys.presX = (2.0*sys.kinX + sys.virX)/md_cell_volume;
  sys.presY = (2.0*sys.kinY + sys.virY)/md_cell_volume;
  sys.presZ = (2.0*sys.kinZ + sys.virZ)/md_cell_volume;	

  /* current pressure [ = Trace[Pab]/3 ] */
  sys.pres = (sys.presX + sys.presY + sys.presZ)/3.0;

}

/*****
*   pressure control by forced scaling method
*   see "manual/manual.IEMD".
*****/
void   control_press(int d_step)
{
  static double oldPx = 0.0, oldPy = 0.0, oldPz = 0.0;
  double averaged_Px, averaged_Py, averaged_Pz;
  double target_Px, target_Py, target_Pz;
  static short count = 0;
  double SX, SY, SZ;  /* scale factor */
  static double adjuster = 1.0;  /* [TPa^(-1)] */
  double dpx, dpy, dpz;
  int i;

  count++;

  oldPx += sys.presX;
  oldPy += sys.presY;
  oldPz += sys.presZ;	
	
  if(sys.step % d_step == 0) {

    averaged_Px = oldPx / ((double)count);
    averaged_Py = oldPy / ((double)count);
    averaged_Pz = oldPz / ((double)count);
    oldPx = 0.0;
    oldPy = 0.0;
    oldPz = 0.0;
    count = 0;

    target_Px = ctl.press_X;
    target_Py = ctl.press_Y;
    target_Pz = ctl.press_Z;

    dpx = (averaged_Px - target_Px)*adjuster; /* [TPa] x [TPa^(-1)] */
    dpy = (averaged_Py - target_Py)*adjuster;
    dpz = (averaged_Pz - target_Pz)*adjuster;

    /* evaluate a cell size scaling factor [SX, SY, SZ]        */
    /* where "sys.pres_scalling_factor" is a constant (=0.5).  */
    SX = 1.0 + atan( dpx )*sys.pres_scalling_factor;
    SY = 1.0 + atan( dpy )*sys.pres_scalling_factor;
    SZ = 1.0 + atan( dpz )*sys.pres_scalling_factor;

    /* scaling */
    sys.Lx *= SX;
    sys.Ly *= SY;
    sys.Lz *= SZ;

    /* position shift after the cell size scaling */
    for(i=0; i<sys.N; i++){
      sys.rx[i] *= SX;
      sys.ry[i] *= SY;
      sys.rz[i] *= SZ;
    }

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
      sys.vx[i] *= scale;
      sys.vy[i] *= scale;
      sys.vz[i] *= scale;
    }

   printf("---------------------------------------------------------------\n");
   printf(">> Temperature scaling\n");

  }
}




