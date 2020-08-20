#define XR_INCLUDED
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <X11/Xlib.h>
#include  <X11/Xutil.h>
#include  "../md.h"
#include  "../prototypes.h"
#include  "xr.h"

void   open_xunit(void)   /* XUNIT */
{
  short i=0, ii, r;
  short ptx, pty, ptz;
  short xyzpts_x, xyzpts_y;
  int x_pos, y_pos;
  int x_shift, y_shift, xyz_shift;
  int xyzbox_x, xyzbox_y, roc_x[10], roc_y[10];
  double x_scale, y_scale, z_scale;
  double unit_x, unit_y, unit_z;
  double x, y, z, box_x, box_y, box_z;
  double cos_phy, sin_phy, cos_the, sin_the;

  cos_phy = cos(xunit_phy);   /* [1] Y-rotation */
  sin_phy = sin(xunit_phy);
  cos_the = cos(xunit_theta); /* [2] X-rotation */
  sin_the = sin(xunit_theta);

  /*- Preparate for the Y/X-axis rotation in x-y-z window -*/
  xyz_shift = 25.0+7.0; /* 50.0 = 200.0/2/2 */

  x_shift = (int)(xyz_shift*(cos_phy+sin_phy-1.7));
  y_shift = (int)(xyz_shift*(sin_phy*sin_the+cos_the-
                             cos_phy*sin_the-1.7));

  /* re-scale all of the atomic positions */
  x_scale = 200.0/sys.Lx;
  y_scale = 200.0/sys.Ly;
  z_scale = 200.0/sys.Lz;

  unit_x = 1.1*sys.Lx/(double)sys.nx;
  unit_y = 1.1*sys.Ly/(double)sys.ny;
  unit_z = 1.1*sys.Lz/(double)sys.nz;

  for(z=0.0;z<2.0;z++)
    for(y=0.0;y<2.0;y++)
      for(x=0.0;x<2.0;x++)
	{
	  box_x = x * unit_x * x_scale +9.0;
	  box_y = y * unit_y * y_scale +9.0;
	  box_z = z * unit_z * z_scale +9.0;

	  xyzbox_x = (int)(box_x*cos_phy+box_z*sin_phy)-x_shift;
	  xyzbox_y = (int)(box_x*sin_phy*sin_the+
			   box_y*cos_the-box_z*cos_phy*sin_the)-y_shift;
	  i++;

	  roc_x[i] = (int)(double)((xyzbox_x+5)*xunit_zoom)+xunit_shift_x;
	  roc_y[i] = (int)(double)((xyzbox_y+5)*xunit_zoom)+xunit_shift_y;
	}  

  XSetForeground(d, gc, color[White]);

  XDrawLine(d, w[6], gc, roc_x[1], roc_y[1], roc_x[2], roc_y[2]);
  XDrawLine(d, w[6], gc, roc_x[3], roc_y[3], roc_x[4], roc_y[4]);
  XDrawLine(d, w[6], gc, roc_x[5], roc_y[5], roc_x[6], roc_y[6]);
  XDrawLine(d, w[6], gc, roc_x[7], roc_y[7], roc_x[8], roc_y[8]);

  XDrawLine(d, w[6], gc, roc_x[1], roc_y[1], roc_x[5], roc_y[5]);
  XDrawLine(d, w[6], gc, roc_x[2], roc_y[2], roc_x[6], roc_y[6]);
  XDrawLine(d, w[6], gc, roc_x[3], roc_y[3], roc_x[7], roc_y[7]);
  XDrawLine(d, w[6], gc, roc_x[4], roc_y[4], roc_x[8], roc_y[8]);

  XDrawLine(d, w[6], gc, roc_x[1], roc_y[1], roc_x[3], roc_y[3]);
  XDrawLine(d, w[6], gc, roc_x[2], roc_y[2], roc_x[4], roc_y[4]);
  XDrawLine(d, w[6], gc, roc_x[5], roc_y[5], roc_x[7], roc_y[7]);
  XDrawLine(d, w[6], gc, roc_x[6], roc_y[6], roc_x[8], roc_y[8]);

  xunit_drag();
  network();  /* network resolver */

  for(i=0; i<sys.N; i++) {
    if(sys.rx[i]<unit_x && sys.ry[i]<unit_y && sys.rz[i]<unit_z){

      ptx = (int)(sys.rx[i]*x_scale)+9;
      pty = (int)(sys.ry[i]*y_scale)+9;
      ptz = (int)(sys.rz[i]*z_scale)+9;

      xyzpts_x = (int)((double)ptx*cos_phy+(double)ptz*sin_phy)-x_shift;
      xyzpts_y = (int)((double)ptx*sin_phy*sin_the+
		     (double)pty*cos_the-(double)ptz*cos_phy*sin_the)-y_shift;
    
      x_pos = (int)(double)((xyzpts_x+5)*xunit_zoom)+xunit_shift_x;
      y_pos = (int)(double)((xyzpts_y+5)*xunit_zoom)+xunit_shift_y;

      for(ii=0;ii<GRADATION;ii++){
	r = GRADATION*2 - 2*ii;
	if(sys.ion[i]<7){
	  XSetForeground(d, gc, Col[sys.ion[i]*7+ii].pixel);
	} 
	else{
	  XSetForeground(d, gc, Col[7*7+ii].pixel);
	}
	XFillArc(d, w[6], gc, x_pos-r, y_pos-r, r*2, r*2,0, 360*64);
      }

      XSetForeground(d, gc, color[White]);
      sprintf(stringbuffer, "%d", i);
      XDrawString(d, w[6], gc, x_pos-5, y_pos-15,
		  stringbuffer, strlen(stringbuffer));
    }
  }
}

void xunit_drag(void)
{
  static int x0=0;
  static int y0=0;
  static int button_pressed=0;

  XSelectInput(d, w[6], ButtonPressMask | Button1MotionMask |
	       Button3MotionMask | ButtonReleaseMask);
  if(XCheckWindowEvent(d, w[6], ButtonPressMask, &e)){
    x0 = e.xbutton.x;  
    y0 = e.xbutton.y; 
    if(e.xbutton.button==1) button_pressed = 1;
    if(e.xbutton.button==2) button_pressed = 2;
    if(e.xbutton.button==3) button_pressed = 3;
    printf("button %d pressed\n",button_pressed);
  }

  if(button_pressed==1){
    XCheckWindowEvent(d, w[6], Button1MotionMask, &e);
    /*    printf("button1 moved\n"); */
    xunit_phy += 2.0*(double)(e.xbutton.x-x0)*PI/400.0;
    xunit_theta += 2.0*(double)(e.xbutton.y-y0)*PI/400.0;
    x0 = e.xbutton.x;  
    y0 = e.xbutton.y; 
  }

  if(button_pressed==3){
    XCheckWindowEvent(d, w[6], Button3MotionMask, &e);
    /*    printf("button3 moved\n"); */
    xunit_shift_x += (e.xbutton.x-x0);
    xunit_shift_y += (e.xbutton.y-y0);
    x0 = e.xbutton.x;  
    y0 = e.xbutton.y; 
  }

  if(XCheckWindowEvent(d, w[6], ButtonReleaseMask, &e)){
    button_pressed = 0;
    /*    printf("released\n"); */
  }
}
