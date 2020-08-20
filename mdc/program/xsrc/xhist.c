#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <X11/Xlib.h>
#include  <X11/Xutil.h>
#include  "../md.h"
#include  "../prototypes.h"
#include  "xr.h"

void open_xhist(void)
{
  int i;
  XPoint Ek, Ep, Et;
  XPoint temp, pres;

  Ek.x = (int)(300*sys.step/ctl.calc_max)+80;
  Ek.y = (int)(-60*(sys.kin*sys.perMol)/(1000*EK_MAX))+140;
  Ep.x = (int)(300*sys.step/ctl.calc_max)+80;
  Ep.y = (int)(60*(sys.pot*sys.perMol)/(1000*EP_MAX))+140;
  Et.x = (int)(300*sys.step/ctl.calc_max)+80;
  Et.y = (int)(60*(sys.kin*sys.e2t+sys.pot*sys.perMol)/(1000*ET_MAX))+200;
  temp.x = (int)(300*sys.step/ctl.calc_max)+80;
  temp.y = (int)(-60*(sys.kin*sys.e2t)/TEMP_MAX)+320;
  pres.x = (int)(300*sys.step/ctl.calc_max)+80;
  pres.y = (int)(-60*(sys.pres*sys.pp2gpa)/PRES_MAX)+380;

  if(Ek.y<80)      Ek.y+=60;
  if(Ek.y>140)     Ek.y-=60;
  if(Ep.y<140)     Ep.y+=60;
  if(Ep.y>200)     Ep.y-=60;
  if(Et.y<200)     Et.y+=60;
  if(Et.y>260)     Et.y-=60;
  if(temp.y<260)   temp.y+=60;
  if(temp.y>320)   temp.y-=60;
  if(pres.y<320)   pres.y+=60;
  if(pres.y>380)   pres.y-=60;

  XSetForeground(d, gc, color[Black]);
  XClearArea(d, w[5], 0, 0, 80, 400, True); 

  for(i=0;i<5;i++){
    XDrawRectangle(d, w[5], gc, 80, 80+60*i, 300, 60); 
  }
  XDrawString (d, w[5], gc, 10,  90, "Ek[kJ/mol]", 10);
  XDrawString (d, w[5], gc, 10, 150, "Ep[kJ/mol]", 10);
  XDrawString (d, w[5], gc, 10, 210, "Et[kJ/mol]", 10);
  XDrawString (d, w[5], gc, 10, 270, "T[K]", 4);
  XDrawString (d, w[5], gc, 10, 330, "P[GPa]", 6);

  /* Time Step */
  sprintf(stringbuffer, "Step = %d", sys.step);
  XDrawString(d, w[5], gc, 10, 20, stringbuffer, strlen(stringbuffer));

  /* Kinetic Energy */
  sprintf(stringbuffer, "%.2f", sys.kin*sys.perMol/1000.0);
  XDrawString(d, w[5], gc, 10, 110, stringbuffer, strlen(stringbuffer));

  /* Potential Energy */
  sprintf(stringbuffer, "%.2f", sys.pot*sys.perMol/1000.0);
  XDrawString(d, w[5], gc, 10, 170, stringbuffer, strlen(stringbuffer));

  /* Total Energy */
  sprintf(stringbuffer, "%.2f", (sys.kin+sys.pot)*sys.perMol/1000.0);
  XDrawString(d, w[5], gc, 10, 230, stringbuffer, strlen(stringbuffer));

  /* Temperature display */
  sprintf(stringbuffer, "%6.2f", sys.kin*sys.e2t);
  XDrawString(d, w[5], gc, 10, 290, stringbuffer, strlen(stringbuffer));

  /* Pressure display */
  sprintf(stringbuffer, "%6.3f", sys.pres*sys.pp2gpa);
  XDrawString(d, w[5], gc, 10, 350, stringbuffer, strlen(stringbuffer));

  XSetForeground(d, gc, color[Red]);
  XDrawPoint(d, w[5], gc, Ek.x, Ek.y);
  XDrawPoint(d, w[5], gc, Ep.x, Ep.y);
  XDrawPoint(d, w[5], gc, Et.x, Et.y);
  XDrawPoint(d, w[5], gc, temp.x, temp.y);
  XDrawPoint(d, w[5], gc, pres.x, pres.y);
  
}
