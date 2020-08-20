#define XR_INCLUDED
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <X11/Xlib.h>
#include  <X11/Xutil.h>
#include  "../md.h"
#include  "../prototypes.h"
#include  "xr.h"

void xd_init_set(void)  /* main display set */
{
  short i, ii;

  d = XOpenDisplay(NULL); /* connect to the server */

  /* Set Color in XD-NET display. See also the "enum" declaration */
  color[Red]         = MyColor(d,"red");
  color[Blue]        = MyColor(d,"blue");
  color[Green]       = MyColor(d,"green");
  color[Yellow]      = MyColor(d,"yellow");
  color[Orange]      = MyColor(d,"orange");
  color[Cyan]        = MyColor(d,"cyan");
  color[BlueViolet]  = MyColor(d,"blue violet");
  color[LimeGreen]   = MyColor(d,"lime green");
  color[Coral]       = MyColor(d,"coral");
  color[Khaki]       = MyColor(d,"khaki");
  color[SpringGreen] = MyColor(d,"spring green");
  color[LightGray]   = MyColor(d,"light gray");
  color[Magenta]     = MyColor(d,"magenta");
  color[SlateBlue]   = MyColor(d,"slate blue");
  color[Black]       = MyColor(d,"black");
  color[White]       = MyColor(d,"White");

  /* Window Position Set */
  w[0] = XCreateSimpleWindow(d, RootWindow(d,0), 0, 0, 480, 500, 2,
			     BlackPixel(d, 0), WhitePixel(d, 0)); /* root */
  w[1] = XCreateSimpleWindow(d, w[0],  30,  40, 200, 200, 2, 
			     BlackPixel(d, 0), WhitePixel(d, 0)); /* x-y  */
  w[2] = XCreateSimpleWindow(d, w[0], 250,  40, 200, 200, 2,
			     BlackPixel(d, 0), WhitePixel(d, 0)); /* y-z  */
  w[3] = XCreateSimpleWindow(d, w[0],  30, 270, 200, 200, 2, 
			     BlackPixel(d, 0), WhitePixel(d, 0)); /* z-x  */
  w[4] = XCreateSimpleWindow(d, w[0], 250, 270, 200, 200, 2,
			       BlackPixel(d, 0), WhitePixel(d, 0)); /* zyx  */

  /* New Windows */
  w[5] = XCreateSimpleWindow(d, RootWindow(d,0), 500, 0, 400, 400, 2,
			     BlackPixel(d, 0), WhitePixel(d, 0)); /* XHIST */
  w[6] = XCreateSimpleWindow(d, RootWindow(d,0), 500, 0, 445, 450, 2,
			     WhitePixel(d, 0), BlackPixel(d, 0)); /* XUNIT */
  w[7] = XCreateSimpleWindow(d, RootWindow(d,0), 500, 0, 600, 600, 2,
			     BlackPixel(d, 0), WhitePixel(d, 0)); /* XNET */
  w[8] = XCreateSimpleWindow(d, RootWindow(d,0), 500, 0, 400, 400, 2,
			     BlackPixel(d, 0), WhitePixel(d, 0)); /* XMSD */

  /* Icon Position Set at XD root Window */
  for(i=0;i<12;i++){
    ic[i] = XCreateSimpleWindow(d, w[0], 40*i, 2, 40, 20, 2,
				BlackPixel(d, 0), WhitePixel(d, 0));
  }

  /* Set Color GRADATION for XUNIT display */
  cmap = DefaultColormap(d,0);

  for(ii=0;ii<GRADATION;ii++){
    Col[ii].red  = 50000;
    Col[ii].blue  = 65535 *ii /GRADATION; 
    Col[ii].green = 65535 *ii /GRADATION;
    XAllocColor(d, cmap, &Col[ii]);
  }
  for(ii=GRADATION;ii<GRADATION*2;ii++){
    Col[ii].red  = 65535 *(ii-GRADATION) /GRADATION;
    Col[ii].blue  = 50000;
    Col[ii].green = 65535 *(ii-GRADATION) /GRADATION;
    XAllocColor(d, cmap, &Col[ii]);
  }
  for(ii=GRADATION*2;ii<GRADATION*3;ii++){
    Col[ii].red  = 65535 *(ii-GRADATION*2) /GRADATION;
    Col[ii].blue  = 65535 *(ii-GRADATION*2) /GRADATION;
    Col[ii].green = 50000;
    XAllocColor(d, cmap, &Col[ii]);
  }
  for(ii=GRADATION*3;ii<GRADATION*4;ii++){
    Col[ii].red  = 65535 *(ii-GRADATION*3) /GRADATION;
    Col[ii].blue  = 50000;
    Col[ii].green  = 50000;
    XAllocColor(d, cmap, &Col[ii]);
  }
  for(ii=GRADATION*4;ii<GRADATION*5;ii++){
    Col[ii].red  = 50000;
    Col[ii].blue  = 65535 *(ii-GRADATION*4) /GRADATION;
    Col[ii].green  = 50000;
    XAllocColor(d, cmap, &Col[ii]);
  }
  for(ii=GRADATION*5;ii<GRADATION*6;ii++){
    Col[ii].red  = 50000;
    Col[ii].blue  = 50000;
    Col[ii].green  = 65535 *(ii-GRADATION*5) /GRADATION;
    XAllocColor(d, cmap, &Col[ii]);
  }
  for(ii=GRADATION*6;ii<GRADATION*7;ii++){
    Col[ii].red  = 65535 *(ii-GRADATION*6) /GRADATION;
    Col[ii].green  = 65535 *(ii-GRADATION*6) /GRADATION;
    Col[ii].blue  = 65535 *(ii-GRADATION*6) /GRADATION;
    XAllocColor(d, cmap, &Col[ii]);
  }

  /* change the Window attributes to a BackingStore type */
  attr.backing_store = Always;
  for(i=0;i<9;i++){
    XChangeWindowAttributes(d, w[i], CWBackingStore, &attr );
  }
  for(i=0;i<12;i++){
    XChangeWindowAttributes(d, ic[i], CWBackingStore, &attr );
    XSelectInput(d, ic[i], ExposureMask|ButtonPressMask|ButtonReleaseMask);
  }

  /* Mapping */
  XMapWindow(d, w[0]);
  XMapSubwindows(d, w[0]);
  gc = XCreateGC (d, w[0], 0, 0);

  /* Make Titles */
  XStoreName(d, w[0], "IEMD");
  XStoreName(d, w[5], "IEMD HISTORY");
  XStoreName(d, w[6], "IEMD UNIT CELL");
  XStoreName(d, w[7], "IEMD NETWORK");
  XStoreName(d, w[8], "IEMD MSD");
  XSetIconName(d, w[0], "IEMD");
  XSetIconName(d, w[5], "XHIST");
  XSetIconName(d, w[6], "XUNIT");
  XSetIconName(d, w[7], "XNET");
  XSetIconName(d, w[8], "XMSD");

  XSetForeground(d, gc, color[Black]);
  XDrawString(d, w[0], gc, 120, 260, "x-y", 3);
  XDrawString(d, w[0], gc, 330, 260, "y-z", 3);
  XDrawString(d, w[0], gc, 120, 490, "z-x", 3);
  XDrawString(d, w[0], gc, 330, 490, "xyz", 3);

  /* Flames of XD window */    
  XDrawRectangle(d, w[0], gc,  30-3,  40-3, 200+9, 200+9);
  XDrawRectangle(d, w[0], gc, 250-3,  40-3, 200+9, 200+9);
  XDrawRectangle(d, w[0], gc,  30-3, 270-3, 200+9, 200+9);
  XDrawRectangle(d, w[0], gc, 250-3, 270-3, 200+9, 200+9);

  /* memory allocation to use in XNET (the structure is defineded in md.h) */
  net.pos = (short *)calloc(10*sys.N+2, sizeof(short)); 
  net.x_pos = (int *)calloc(sys.N+2, sizeof(int)); 
  net.y_pos = (int *)calloc(sys.N+2, sizeof(int)); 

  /* XNET parameter initialization */
  xnet_theta=0.0;    /* X axis rotation [radian] */
  xnet_phy=0.0;      /* Y axis rotation [radian] */
  xnet_zoom=1.5;     /* ZOOM controller */
  xnet_shift_x=30;    /* X axis shift [pixel] */
  xnet_shift_y=30;    /* Y axis shift [pixel] */

  /* XUNIT parameter initialization */
  xunit_theta=0.0;    /* X axis rotation [radian] */
  xunit_phy=0.0;      /* Y axis rotation [radian] */
  xunit_zoom=4.0;     /* ZOOM controller */
  xunit_shift_x=-60;    /* X axis shift [pixel] */
  xunit_shift_y=-70;    /* Y axis shift [pixel] */
}



