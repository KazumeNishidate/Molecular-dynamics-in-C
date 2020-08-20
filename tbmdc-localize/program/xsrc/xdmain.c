#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "../md.h"
#include  "../prototypes.h"
#include  <X11/Xlib.h>
#include  <X11/Xutil.h>
#include  "xr.h"

void xd(void)
{
  static short xd_cnt=0;

  if(xd_cnt==0){ /* initialize the XD window system */
    xd_init_set();
    xd_cnt = 1;
  }

  open_xd();     /* scale and plot : XD */

  xd_icon_set(); /* check the icon event : XD XUNIT XNET XHIST XMSD */

  open_xhist();  /* store the MD history data to display in XHIST */

  open_xmsd();   /* store the MSD data to display in XMSD */

  XFlush(d);
}

void open_xd(void)
{
  short  i, ptx, pty, ptz;
  double x_scale, y_scale, z_scale;

  static short counter=0;
  static int x_shift, y_shift, xyz_shift;
  static double xd_theta=1.3, xd_phy=1.3;  /* [radian] */
  static double xd_cos_phy, xd_sin_phy, xd_cos_the, xd_sin_the;

  XPoint xypts[XD_MAXPTS], yzpts[XD_MAXPTS];
  XPoint zxpts[XD_MAXPTS], xyzpts[XD_MAXPTS];

  /*- preparation for the Y/X-axis rotation in x-y-z window -*/
  if(counter == 0){
    /* xd_theta: X rotation angle in [radian] for XYZ */
    /* xd_phy:   Y rotation angle in [radian] for XYZ */

    xd_cos_phy = cos(xd_phy);   /* [1] Y-rotation */
    xd_sin_phy = sin(xd_phy);
    xd_cos_the = cos(xd_theta); /* [2] X-rotation */
    xd_sin_the = sin(xd_theta);

    xyz_shift = 183.0/2.0+7.0;
    x_shift = ((int)xyz_shift*(int)(xd_cos_phy+xd_sin_phy-1.0))+15;
    y_shift = ((int)xyz_shift*(int)(xd_sin_phy*xd_sin_the+xd_cos_the-
			       xd_cos_phy*xd_sin_the-1.0))-10;
    counter = 1;
  }

  /*- re-scale all of the atomic positions -*/
  x_scale = (int)(180.0/sys.Lx);
  y_scale = (int)(180.0/sys.Ly);
  z_scale = (int)(180.0/sys.Lz);      

  /*- XD projection views for [XY YZ ZX XYZ] surfaces -*/
  for(i=0; i<sys.N; i++) {
    ptx = (int)(sys.rx[i]*x_scale)+9; 
    pty = (int)(sys.ry[i]*y_scale)+9; 
    ptz = (int)(sys.rz[i]*z_scale)+9; 
    
    xypts[i].x = ptx;
    xypts[i].y = pty;
    yzpts[i].x = pty;
    yzpts[i].y = ptz;
    zxpts[i].x = ptz;
    zxpts[i].y = ptx;
    xyzpts[i].x = (int)((double)ptx*xd_cos_phy+
			(double)ptz*xd_sin_phy)-x_shift;
    xyzpts[i].y = (int)((double)ptx*xd_sin_phy*xd_sin_the+
		(double)pty*xd_cos_the-(double)ptz*xd_cos_phy*xd_sin_the)-
		  y_shift;

    XSetForeground (d, gc, color[sys.ion[i]]);
    XDrawPoint(d, w[1], gc, xypts[i].x, xypts[i].y);
    XDrawPoint(d, w[2], gc, yzpts[i].x, yzpts[i].y);
    XDrawPoint(d, w[3], gc, zxpts[i].x, zxpts[i].y);
    XDrawPoint(d, w[4], gc, xyzpts[i].x, xyzpts[i].y);
  }
}
/********************************************************************/
 
