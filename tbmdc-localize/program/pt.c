#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*
*  calculation of pressure tensor P_xx, P_yy, and P_zz using "virial".
*/

void   calc_press(void)
{
  double md_cell_volume;

  md_cell_volume = sys.Lx*sys.Ly*sys.Lz;

  sys.presX = (2.0*sys.kinX + press.virX/6.0)/md_cell_volume;
  sys.presY = (2.0*sys.kinY + press.virY/6.0)/md_cell_volume;
  sys.presZ = (2.0*sys.kinZ + press.virZ/6.0)/md_cell_volume;	

  /* current pressure     */
  sys.pres = (sys.presX + sys.presY + sys.presZ)/3.0;
}

/*****
*   pressure control by forced scaling method
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

/* evaluate a cell size scaling factor [SX, SY, SZ]        */
/* where "sys.pres_scalling_factor" is a constant (=0.5).  */
    SX = 1.0 + atan(averaged_Px - target_Px)*sys.pres_scalling_factor;
    SY = 1.0 + atan(averaged_Py - target_Py)*sys.pres_scalling_factor;
    SZ = 1.0 + atan(averaged_Pz - target_Pz)*sys.pres_scalling_factor;

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
/*
   printf("---------------------------------------------------------------\n");
   printf(">> Pressure scaling\n");
*/
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

  if(sys.step % d_step == 0) {
    averaged_K = oldk / (double) count;
    oldk = 0.0;
    count = 0;
    d_step=0;
  }

  if(d_step == 0) {
    tempK  = temp/sys.e2t;              /* [energy] -> temperature   */
    scale = sqrt(tempK / averaged_K);   /* evaluate a scaling factor */
    for(i=0; i< sys.N; i++) {
      sys.vx[i] *= scale;
      sys.vy[i] *= scale;
      sys.vz[i] *= scale;
    }
/*
   printf("---------------------------------------------------------------\n");
   printf(">> Temperature scaling\n");
*/
  }
}




