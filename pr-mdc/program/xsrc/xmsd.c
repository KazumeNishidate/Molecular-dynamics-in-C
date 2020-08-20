#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <X11/Xlib.h>
#include  <X11/Xutil.h>
#include  "../md.h"
#include  "../prototypes.h"
#include  "xr.h"

void open_xmsd(void)
{
  int i, msd_max=10;
  static int xmsd_cnt=0;
  static XPoint msd_ion[1];

  if(xmsd_cnt==0){
    XSetForeground(d, gc, color[Black]);
    XDrawRectangle(d, w[8], gc, 20, 20, 360, 360);
    XDrawLine(d, w[8], gc, 20, 200, 380, 200);

    xmsd_cnt = 1;
  }

  XClearArea(d, w[8], 60, 30, 50, 40+20*(ctl.kinds_of_ions-2), True); 
  XClearArea(d, w[8], 0,0, 50, 15, True); 

  for(i=0;i<ctl.kinds_of_ions;i++){

    msd_ion[0].x = (int)(360*sys.step/ctl.calc_max)+20;
    msd_ion[0].y = -(int)(360*msd.value[i]/10)+380;

    if(msd_ion[0].y<20){
      msd_ion[0].y = -(int)(360*msd.value[i]/10)+380+360;
      msd_max+=10;
    }
    if(msd_ion[0].y>380){
      msd_ion[0].y = -(int)(360*msd.value[i]/10)+380-360;
      msd_max-=10;
    }

    XSetForeground(d, gc, color[i]);
    XDrawPoint(d, w[8], gc, msd_ion[0].x, msd_ion[0].y);

    sprintf(stringbuffer, "ion%d=%.3f",i,msd.value[i]);
    XDrawString(d, w[8], gc, 30, 40+20*i, stringbuffer, strlen(stringbuffer));
  }    
  XSetForeground(d, gc, color[Black]);

  sprintf(stringbuffer, "MAX=%d",msd_max);
  XDrawString(d, w[8], gc, 5, 15, stringbuffer, strlen(stringbuffer));

}
