#define XR_INCLUDED
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <X11/Xlib.h>
#include  <X11/Xutil.h>
#include  "../md.h"
#include  "../prototypes.h"
#include  "xr.h"

void   open_xnet(void)   /* XNET */
{
  short i, ii, jj;
  short ptx, pty, ptz;
  double x_scale, y_scale, z_scale;
  int x_shift, y_shift, xyz_shift;
  short xyzpts_x, xyzpts_y;
  double cos_phy, sin_phy, cos_the, sin_the;

  xnet_drag(); /* interactively change the net.length value using Mouse :-) */
  network();   /* network resolver */

  cos_phy = cos(xnet_phy);   /* [1] Y-rotation */
  sin_phy = sin(xnet_phy);
  cos_the = cos(xnet_theta); /* [2] X-rotation */
  sin_the = sin(xnet_theta);

  /*- Preparate for the Y/X-axis rotation in x-y-z window -*/
  xyz_shift = 100.0+7.0; /* 100.0 = 200.0/2 */

  x_shift = (int)(xyz_shift*(cos_phy+sin_phy-1.7));
  y_shift = (int)(xyz_shift*(sin_phy*sin_the+cos_the-
                             cos_phy*sin_the-1.7));

  /* re-scale all of the atomic positions */
  x_scale = 200.0/sys.Lx;
  y_scale = 200.0/sys.Ly;
  z_scale = 200.0/sys.Lz;

  for(i=0; i<sys.N; i++) {
    ptx = (int)(sys.rx[i]*x_scale)+9;
    pty = (int)(sys.ry[i]*y_scale)+9;
    ptz = (int)(sys.rz[i]*z_scale)+9;
    
    xyzpts_x = (int)((double)ptx*cos_phy+
		     (double)ptz*sin_phy)-x_shift;
    xyzpts_y = (int)((double)ptx*sin_phy*sin_the+
		     (double)pty*cos_the-
		     (double)ptz*cos_phy*sin_the)-y_shift;

    net.x_pos[i] = (int)((double)(xyzpts_x)*xnet_zoom)+xnet_shift_x;
    net.y_pos[i] = (int)((double)(xyzpts_y)*xnet_zoom)+xnet_shift_y;

    XSetForeground(d, gc, color[sys.ion[i]]);
    XFillArc(d, w[7], gc, net.x_pos[i]-5, net.y_pos[i]-5, 10, 10, 0, 360*64);
  }

  ii=0;
  XSetForeground(d, gc, color[Black]);
  for(jj=0;jj<net.counter/2;jj++){
    XDrawLine(d, w[7], gc,
	      net.x_pos[net.pos[ii]], net.y_pos[net.pos[ii]],
	      net.x_pos[net.pos[ii+1]], net.y_pos[net.pos[ii+1]]);
    ii+=2;
  }
}

/***
*  Interactively change the net.length value using Mouse button.
*
*  The "net.length" is defined in network_resolver structure of 
*  "program/md.h".  see also "network.c".
***/
void xnet_drag(void)   
{
  static int x0=0;
  static int y0=0;
  static int button_pressed=0; 

  XSelectInput(d, w[7], ButtonPressMask | ButtonReleaseMask);

  if(XCheckWindowEvent(d, w[7], ButtonPressMask, &e)){
    x0 = e.xbutton.x;  
    y0 = e.xbutton.y; 
    if(e.xbutton.button==1) button_pressed = 1;
    if(e.xbutton.button==2) button_pressed = 2;
    if(e.xbutton.button==3) button_pressed = 3;
    /* printf("button %d pressed\n",button_pressed); */
  }

  if(button_pressed==1){
    net.length+=0.5;     /* increase the net.length by 0.5 [A] */
  }
  if(button_pressed==3){
    net.length-=0.5;     /* decrease the net.length by 0.5 [A] */
  }
  if(button_pressed==2){
    net.length=10.0;     /* reset the network length */
  }
  if(XCheckWindowEvent(d, w[7], ButtonReleaseMask, &e)){
    button_pressed = 0;  /* reset the button FLAG */
  }

}

