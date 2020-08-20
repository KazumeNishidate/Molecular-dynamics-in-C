#define XR_INCLUDED
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <X11/Xlib.h>
#include  <X11/Xutil.h>
#include  "../md.h"
#include  "../prototypes.h"
#include  "xr.h"

void xd_icon_set(void)
{
  static int counter = 0;

  if(counter==0){ /* Icon title set */
    XSetForeground (d, gc, color[Black]);
    XDrawString(d, ic[0],  gc, 8, 15, "QUIT",  4);  
    XDrawString(d, ic[1],  gc, 5, 15, "PAUSE", 5);  
    XDrawString(d, ic[2],  gc, 5, 15, "CLEAR", 5);  
    XDrawString(d, ic[3],  gc, 8, 15, "HEAT",  4);  
    XDrawString(d, ic[4],  gc, 8, 15, "COOL",  4);  
    XDrawString(d, ic[5],  gc, 5, 15, "XHIST", 5);  
    XDrawString(d, ic[6],  gc, 5, 15, "XUNIT", 5);  
    XDrawString(d, ic[7],  gc, 8, 15, "XNET",  4);  
    XDrawString(d, ic[8],  gc, 8, 15, "XMSD",  4);  
    XDrawString(d, ic[9],  gc, 5, 15, "-----", 5);  
    XDrawString(d, ic[10], gc, 5, 15, "-----", 5);  
    XDrawString(d, ic[11], gc, 5, 15, "-----", 5);  
  }

  /* Icon's Function Set */
  /*-------------------- QUIT BUTTON  (ic[0]) -------------------------*/
  if(XCheckWindowEvent(d, ic[0], ButtonPressMask, &e)){
    exit(0);
  }

  /*-------------------- PAUSE BUTTON  (ic[1]) ------------------------*/
  if(XCheckWindowEvent(d, ic[1], ButtonPressMask, &e)){
    XMoveResizeWindow(d, ic[1], 40, 0, 40, 20);
    while(1){
      if(XCheckWindowEvent(d, ic[0], ButtonPressMask, &e)) exit(0);
      if(XCheckWindowEvent(d, ic[1], ButtonPressMask,&e))	{
	XMoveResizeWindow(d, ic[1], 40, 2, 40, 20);	
	break;
      }
      xd_icon_clear();
    }
  }
  XClearArea(d,w[5],0,0,80,400,True);  
  XClearWindow(d,w[6]); /* for xroll always display */
  XClearWindow(d,w[7]); /* for xnet always display */

  xd_icon_clear();  /* check the CLEAR button     */
  xd_icon_heat();   /* check the HEAT up button   */
  xd_icon_cool();   /* check the COOL down button */
  xd_icon_xhist();  /* check the XHIST button     */
  xd_icon_xunit();  /* check the XUNIT button     */
  xd_icon_xnet();   /* check the XNET button      */
  xd_icon_xmsd();   /* check the XMSD button      */
}
/*-------------------- CLEAR BUTTON  (ic[2]) ------------------------*/

void xd_icon_clear(void)  
{
  if(XCheckWindowEvent(d, ic[2], ButtonPressMask, &e)){
    XSetForeground(d, gc, color[White]);
    XFillRectangle(d, w[1], gc, 0, 0, 200, 200);
    XFillRectangle(d, w[2], gc, 0, 0, 200, 200);
    XFillRectangle(d, w[3], gc, 0, 0, 200, 200);
    XFillRectangle(d, w[4], gc, 0, 0, 200, 200);
  }
}

/*-------------------- HEAT BUTTON  (ic[3]) --------------------------*/
void xd_icon_heat(void)  
{
  static short temp_heat=0;
  switch(temp_heat){
  case 0:
    if(XCheckWindowEvent(d, ic[3], ButtonPressMask, &e)){
      temp_heat = 1;
      XMoveResizeWindow(d, ic[3], 120, 0, 40, 20);
    }
    break;
  case 1:
    ctl.temp += 5;
    if(XCheckWindowEvent(d, ic[3], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic[3], 120, 2, 40, 20);
      temp_heat = 0;
    }
    break;
  }
}  

/*------------------- COOL BUTTON  (ic[4]) ---------------------------*/
void xd_icon_cool(void)  
{
  static short temp_cool=0;

  switch(temp_cool){
  case 0:
    if(XCheckWindowEvent(d, ic[4], ButtonPressMask, &e)){
      temp_cool = 1;
      XMoveResizeWindow(d, ic[4], 160, 0, 40, 20);
    }
    break;
  case 1:
    ctl.temp -= 5;
    if(XCheckWindowEvent(d, ic[4], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic[4], 160, 2, 40, 20);
      temp_cool = 0;
    }
    break;
  }
}

/*------------------ XHIST BUTTON  (ic[5]) ---------------------------*/
void xd_icon_xhist(void)
{
  static short icon_xhist=0;

  switch(icon_xhist){
  case 0:
    if(XCheckWindowEvent(d, ic[5], ButtonPressMask, &e)){
      icon_xhist = 1;
      XMoveResizeWindow(d, ic[5], 200, 0, 40, 20);
      XMapWindow(d, w[5]);
    }
    break;
  case 1:
    open_xhist();
    if(XCheckWindowEvent(d, ic[5], ButtonPressMask, &e)){
      icon_xhist = 0;
      XMoveResizeWindow(d, ic[5], 200, 2, 40, 20);
      XUnmapWindow(d, w[5]);
    }
    break;
  }
}

/*------------------ XUNIT BUTTON  (ic[6]) ---------------------------*/
void xd_icon_xunit(void)
{
  int i;
  static short icon_xunit=0;
  static short icon_xunit_ic0=0, icon_xunit_ic1=0; 
  static short icon_xunit_ic2=0, icon_xunit_ic3=0; 
  static short icon_xunit_ic4=0, icon_xunit_ic5=0; 
  static short icon_xunit_ic6=0, icon_xunit_ic7=0; 
  static short icon_xunit_ic8=0, icon_xunit_ic9=0; 

  switch(icon_xunit){
  case 0:
    if(XCheckWindowEvent(d, ic[6], ButtonPressMask, &e)){
      icon_xunit = 1;
      XMoveResizeWindow(d, ic[6], 240, 0, 40, 20);
      /*** XUNIT ICON open check     ***/

      XSetForeground(d, gc, color[Black]);      
      for(i=0;i<XD_XUNIT_ICONS;i++){
	ic_xunit[i] = XCreateSimpleWindow(d, w[6], 40*i, 2, 40, 20, 2,
					 BlackPixel(d, 0), WhitePixel(d, 0));
      }
      /* XD_NET_ICONS  -->  XD_XUNIT_ICONS   */
      for(i=0;i<XD_XUNIT_ICONS;i++){
	XChangeWindowAttributes(d, ic_xunit[i], CWBackingStore, &attr );
	XSelectInput(d, ic_xunit[i], ExposureMask|
		     ButtonPressMask|ButtonReleaseMask);
      }
      XDrawString(d, ic_xunit[0], gc, 5, 15, "Xrot+", 5);  
      XDrawString(d, ic_xunit[1], gc, 5, 15, "Xrot-", 5);  
      XDrawString(d, ic_xunit[2], gc, 5, 15, "Yrot+", 5);  
      XDrawString(d, ic_xunit[3], gc, 5, 15, "Yrot-", 5);  
      XDrawString(d, ic_xunit[4], gc, 5, 15, "zoom+", 5);  
      XDrawString(d, ic_xunit[5], gc, 5, 15, "zoom-", 5);  
      XDrawString(d, ic_xunit[6], gc, 5, 15, " <== ", 5);  
      XDrawString(d, ic_xunit[7], gc, 5, 15, " ==> ", 5);  
      XDrawString(d, ic_xunit[8], gc, 5, 15, " UP  ", 5);  
      XDrawString(d, ic_xunit[9], gc, 5, 15, "DOWN ", 5);  
      XDrawString(d, ic_xunit[10], gc, 5, 15, "RESET", 5);  
      /*** XUNIT ICON open check end ***/

      XMapSubwindows(d, w[6]);
      XMapWindow(d, w[6]);
      open_xunit();
    }
    break;
  case 1:
    open_xunit();

    /*--------------------------------------------------------*/
    if(icon_xunit_ic0==1 && 
       XCheckWindowEvent(d, ic_xunit[0], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[0], 0, 2, 40, 20);
      icon_xunit_ic0=0;
    }
    if(icon_xunit_ic1==1 && 
       XCheckWindowEvent(d, ic_xunit[1], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[1], 40, 2, 40, 20);
      icon_xunit_ic1=0;
    }
    if(icon_xunit_ic2==1 && 
       XCheckWindowEvent(d, ic_xunit[2], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[2], 80, 2, 40, 20);
      icon_xunit_ic2=0;
    }
    if(icon_xunit_ic3==1 && 
       XCheckWindowEvent(d, ic_xunit[3], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[3], 120, 2, 40, 20);
      icon_xunit_ic3=0;
    }
    if(icon_xunit_ic4==1 && 
       XCheckWindowEvent(d, ic_xunit[4], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[4], 160, 2, 40, 20);
      icon_xunit_ic4=0;
    }
    if(icon_xunit_ic5==1 && 
       XCheckWindowEvent(d, ic_xunit[5], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[5], 200, 2, 40, 20);
      icon_xunit_ic5=0;
    }
    /* check */
    if(icon_xunit_ic6==1 && 
       XCheckWindowEvent(d, ic_xunit[6], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[6], 240, 2, 40, 20);
      icon_xunit_ic6=0;
    }
    if(icon_xunit_ic7==1 && 
       XCheckWindowEvent(d, ic_xunit[7], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[7], 280, 2, 40, 20);
      icon_xunit_ic7=0;
    }
    if(icon_xunit_ic8==1 && 
       XCheckWindowEvent(d, ic_xunit[8], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[8], 320, 2, 40, 20);
      icon_xunit_ic8=0;
    }
    if(icon_xunit_ic9==1 && 
       XCheckWindowEvent(d, ic_xunit[9], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xunit[9], 360, 2, 40, 20);
      icon_xunit_ic9=0;
    }
    /*--------------------------------------------------------*/
    if(XCheckWindowEvent(d, ic_xunit[0], ButtonPressMask, &e)){
      icon_xunit_ic0=1;
      XMoveResizeWindow(d, ic_xunit[0], 0, 0, 40, 20);
    }
    if(icon_xunit_ic0==1){
      xunit_theta+=0.03;
      if(xunit_theta > PI2 ) xunit_theta = 0.0; 
    }
    if(XCheckWindowEvent(d, ic_xunit[1], ButtonPressMask, &e)){
      icon_xunit_ic1=1;
      XMoveResizeWindow(d, ic_xunit[1], 40, 0, 40, 20);
    }
    if(icon_xunit_ic1==1){
      xunit_theta-=0.03;
      if(xunit_theta > PI2 ) xunit_theta = 0.0; 
    }
    if(XCheckWindowEvent(d, ic_xunit[2], ButtonPressMask, &e)){
      icon_xunit_ic2=1;
      XMoveResizeWindow(d, ic_xunit[2], 80, 0, 40, 20);
    }
    if(icon_xunit_ic2==1){
      xunit_phy+=0.03;
      if(xunit_phy > PI2 ) xunit_phy = 0.0; 
    }
    if(XCheckWindowEvent(d, ic_xunit[3], ButtonPressMask, &e)){
      icon_xunit_ic3=1;
      XMoveResizeWindow(d, ic_xunit[3], 120, 0, 40, 20);
    }
    if(icon_xunit_ic3==1){
      xunit_phy-=0.03;
      if(xunit_phy > PI2 ) xunit_phy = 0.0; 
    }
    if(XCheckWindowEvent(d,ic_xunit[4],ButtonPressMask,&e)){
      icon_xunit_ic4=1;
      XMoveResizeWindow(d,ic_xunit[4],160,0,40,20);
    }
    if(icon_xunit_ic4==1){
      xunit_zoom+=0.1;
    }
    if(XCheckWindowEvent(d,ic_xunit[5],ButtonPressMask,&e)){
      icon_xunit_ic5=1;
      XMoveResizeWindow(d,ic_xunit[5],200,0,40,20);
    }
    if(icon_xunit_ic5==1){
      xunit_zoom-=0.1;
    }
    /********************************/
    if(XCheckWindowEvent(d,ic_xunit[6],ButtonPressMask,&e)){
      icon_xunit_ic6=1;
      XMoveResizeWindow(d,ic_xunit[6],240,0,40,20);
    }
    if(icon_xunit_ic6==1){
      xunit_shift_x -= 5;
    }
    if(XCheckWindowEvent(d,ic_xunit[7],ButtonPressMask,&e)){
      icon_xunit_ic7=1;
      XMoveResizeWindow(d,ic_xunit[7],280,0,40,20);
    }
    if(icon_xunit_ic7==1){
      xunit_shift_x += 5;
    }
    if(XCheckWindowEvent(d,ic_xunit[8],ButtonPressMask,&e)){
      icon_xunit_ic8=1;
      XMoveResizeWindow(d,ic_xunit[8],320,0,40,20);
    }
    if(icon_xunit_ic8==1){
      xunit_shift_y -= 5;
    }
    if(XCheckWindowEvent(d,ic_xunit[9],ButtonPressMask,&e)){
      icon_xunit_ic9=1;
      XMoveResizeWindow(d,ic_xunit[9],360,0,40,20);
    }
    if(icon_xunit_ic9==1){
      xunit_shift_y += 5;
    }

    /*******************************/
    /* check */
    if(XCheckWindowEvent(d,ic_xunit[10],ButtonPressMask,&e)){
      /* XUNIT parameter initialization */
      xunit_theta=0.0;    /* X axis rotation [radian] */
      xunit_phy=0.0;      /* Y axis rotation [radian] */
      xunit_zoom=4.0;     /* ZOOM controller */
      xunit_shift_x=-60;    /* X axis shift [pixel] */
      xunit_shift_y=-70;    /* Y axis shift [pixel] */
    }
    /*--------------------------------------------------------*/

    /*================================*/
    if(XCheckWindowEvent(d, ic[6], ButtonPressMask, &e)){
      icon_xunit = 0;
      XMoveResizeWindow(d, ic[6], 240, 2, 40, 20);
      XUnmapWindow(d, w[6]);
    }
    break;
  }
}

/*------------------ XNET BUTTON  (ic[7]) ---------------------------*/
void xd_icon_xnet(void)
{
  int i;
  static short icon_xnet=0; 
  static short icon_xnet_ic0=0, icon_xnet_ic1=0; 
  static short icon_xnet_ic2=0, icon_xnet_ic3=0; 
  static short icon_xnet_ic4=0, icon_xnet_ic5=0; 
  static short icon_xnet_ic6=0, icon_xnet_ic7=0; 
  static short icon_xnet_ic8=0, icon_xnet_ic9=0; 

  switch(icon_xnet){
  case 0:
    if(XCheckWindowEvent(d, ic[7], ButtonPressMask, &e)){
      icon_xnet = 1;
      XMoveResizeWindow(d, ic[7], 280, 0, 40, 20);
      XSetForeground(d, gc, color[Black]);      
      for(i=0;i<XD_XNET_ICONS;i++){
	ic_xnet[i] = XCreateSimpleWindow(d, w[7], 40*i, 2, 40, 20, 2,
					 BlackPixel(d, 0), WhitePixel(d, 0));
      }
      for(i=0;i<XD_XNET_ICONS;i++){
	XChangeWindowAttributes(d, ic_xnet[i], CWBackingStore, &attr );
	XSelectInput(d, ic_xnet[i], ExposureMask|
		     ButtonPressMask|ButtonReleaseMask);
      }
      XDrawString(d, ic_xnet[0], gc, 5, 15, "Xrot+", 5);  
      XDrawString(d, ic_xnet[1], gc, 5, 15, "Xrot-", 5);  
      XDrawString(d, ic_xnet[2], gc, 5, 15, "Yrot+", 5);  
      XDrawString(d, ic_xnet[3], gc, 5, 15, "Yrot-", 5);  
      XDrawString(d, ic_xnet[4], gc, 5, 15, "zoom+", 5);  
      XDrawString(d, ic_xnet[5], gc, 5, 15, "zoom-", 5);  
      XDrawString(d, ic_xnet[6], gc, 5, 15, " <== ", 5);  
      XDrawString(d, ic_xnet[7], gc, 5, 15, " ==> ", 5);  
      XDrawString(d, ic_xnet[8], gc, 5, 15, " UP  ", 5);  
      XDrawString(d, ic_xnet[9], gc, 5, 15, "DOWN ", 5);  
      XDrawString(d, ic_xnet[10], gc, 5, 15, "RESET", 5);  
      XDrawString(d, ic_xnet[11], gc, 5, 15, "CLOSE", 5);

      XMapWindow(d, w[7]);
      XMapSubwindows(d, w[7]);
      open_xnet();
    }
    break;
  case 1:
    open_xnet();
    /*--------------------------------------------------------*/
    if(icon_xnet_ic0==1 && 
       XCheckWindowEvent(d, ic_xnet[0], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[0], 0, 2, 40, 20);
      icon_xnet_ic0=0;
    }
    if(icon_xnet_ic1==1 && 
       XCheckWindowEvent(d, ic_xnet[1], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[1], 40, 2, 40, 20);
      icon_xnet_ic1=0;
    }
    if(icon_xnet_ic2==1 && 
       XCheckWindowEvent(d, ic_xnet[2], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[2], 80, 2, 40, 20);
      icon_xnet_ic2=0;
    }
    if(icon_xnet_ic3==1 && 
       XCheckWindowEvent(d, ic_xnet[3], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[3], 120, 2, 40, 20);
      icon_xnet_ic3=0;
    }
    if(icon_xnet_ic4==1 && 
       XCheckWindowEvent(d, ic_xnet[4], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[4], 160, 2, 40, 20);
      icon_xnet_ic4=0;
    }
    if(icon_xnet_ic5==1 && 
       XCheckWindowEvent(d, ic_xnet[5], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[5], 200, 2, 40, 20);
      icon_xnet_ic5=0;
    }
    if(icon_xnet_ic6==1 && 
       XCheckWindowEvent(d, ic_xnet[6], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[6], 240, 2, 40, 20);
      icon_xnet_ic6=0;
    }
    if(icon_xnet_ic7==1 && 
       XCheckWindowEvent(d, ic_xnet[7], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[7], 280, 2, 40, 20);
      icon_xnet_ic7=0;
    }
    if(icon_xnet_ic8==1 && 
       XCheckWindowEvent(d, ic_xnet[8], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[8], 320, 2, 40, 20);
      icon_xnet_ic8=0;
    }
    if(icon_xnet_ic9==1 && 
       XCheckWindowEvent(d, ic_xnet[9], ButtonPressMask, &e)){
      XMoveResizeWindow(d, ic_xnet[9], 360, 2, 40, 20);
      icon_xnet_ic9=0;
    }
    if(XCheckWindowEvent(d,ic_xnet[11],ButtonPressMask,&e)){
      XMoveResizeWindow(d,ic[7],280,2,40,20);
      XUnmapWindow(d,w[7]);
      icon_xnet = 0;
    }

    /*--------------------------------------------------------*/
    if(XCheckWindowEvent(d, ic_xnet[0], ButtonPressMask, &e)){
      icon_xnet_ic0=1;
      XMoveResizeWindow(d, ic_xnet[0], 0, 0, 40, 20);
    }
    if(icon_xnet_ic0==1){
      xnet_theta+=0.03;
      if(xnet_theta > PI2 ) xnet_theta = 0.0; 
    }
    if(XCheckWindowEvent(d, ic_xnet[1], ButtonPressMask, &e)){
      icon_xnet_ic1=1;
      XMoveResizeWindow(d, ic_xnet[1], 40, 0, 40, 20);
    }
    if(icon_xnet_ic1==1){
      xnet_theta-=0.03;
      if(xnet_theta > PI2 ) xnet_theta = 0.0; 
    }
    if(XCheckWindowEvent(d, ic_xnet[2], ButtonPressMask, &e)){
      icon_xnet_ic2=1;
      XMoveResizeWindow(d, ic_xnet[2], 80, 0, 40, 20);
    }
    if(icon_xnet_ic2==1){
      xnet_phy+=0.03;
      if(xnet_phy > PI2 ) xnet_phy = 0.0; 
    }
    if(XCheckWindowEvent(d, ic_xnet[3], ButtonPressMask, &e)){
      icon_xnet_ic3=1;
      XMoveResizeWindow(d, ic_xnet[3], 120, 0, 40, 20);
    }
    if(icon_xnet_ic3==1){
      xnet_phy-=0.03;
      if(xnet_phy > PI2 ) xnet_phy = 0.0; 
    }
    if(XCheckWindowEvent(d,ic_xnet[4],ButtonPressMask,&e)){
      icon_xnet_ic4=1;
      XMoveResizeWindow(d,ic_xnet[4],160,0,40,20);
    }
    if(icon_xnet_ic4==1){
      xnet_zoom+=0.05;
    }
    if(XCheckWindowEvent(d,ic_xnet[5],ButtonPressMask,&e)){
      icon_xnet_ic5=1;
      XMoveResizeWindow(d,ic_xnet[5],200,0,40,20);
    }
    if(icon_xnet_ic5==1){
      xnet_zoom-=0.05;
    }
    /********************************/
    if(XCheckWindowEvent(d,ic_xnet[6],ButtonPressMask,&e)){
      icon_xnet_ic6=1;
      XMoveResizeWindow(d,ic_xnet[6],240,0,40,20);
    }
    if(icon_xnet_ic6==1){
      xnet_shift_x -= 5;
    }
    if(XCheckWindowEvent(d,ic_xnet[7],ButtonPressMask,&e)){
      icon_xnet_ic7=1;
      XMoveResizeWindow(d,ic_xnet[7],280,0,40,20);
    }
    if(icon_xnet_ic7==1){
      xnet_shift_x += 5;
    }
    if(XCheckWindowEvent(d,ic_xnet[8],ButtonPressMask,&e)){
      icon_xnet_ic8=1;
      XMoveResizeWindow(d,ic_xnet[8],320,0,40,20);
    }
    if(icon_xnet_ic8==1){
      xnet_shift_y -= 5;
    }
    if(XCheckWindowEvent(d,ic_xnet[9],ButtonPressMask,&e)){
      icon_xnet_ic9=1;
      XMoveResizeWindow(d,ic_xnet[9],360,0,40,20);
    }
    if(icon_xnet_ic9==1){
      xnet_shift_y += 5;
    }

    /*******************************/
    if(XCheckWindowEvent(d,ic_xnet[10],ButtonPressMask,&e)){
      /* XNET parameter initialization */
      xnet_theta=0.0;    /* X axis rotation [radian] */
      xnet_phy=0.0;      /* Y axis rotation [radian] */
      xnet_zoom=1.5;     /* ZOOM controller */
      xnet_shift_x=30;    /* X axis shift [pixel] */
      xnet_shift_y=30;    /* Y axis shift [pixel] */
    }
    /*--------------------------------------------------------*/

    if(XCheckWindowEvent(d, ic[7], ButtonPressMask, &e)){
      icon_xnet = 0;
      xnet_phy = 0.0;
      xnet_theta = 0.0;
      XMoveResizeWindow(d, ic[7], 280, 2, 40, 20);
      XUnmapWindow(d, w[7]);
    }
    break;
  }
}
/*-------------------- XMSD BUTTON  (ic[8]) --------------------------*/

void xd_icon_xmsd(void)
{
  static short icon_xmsd=0;

  switch(icon_xmsd){
  case 0:
    if(XCheckWindowEvent(d, ic[8], ButtonPressMask, &e)){
      icon_xmsd = 1;
      XMoveResizeWindow(d, ic[8], 320, 0, 40, 20);
      XMapWindow(d, w[8]);
    }
    break;
  case 1:
    open_xmsd();
    if(XCheckWindowEvent(d, ic[8], ButtonPressMask, &e)){
      icon_xmsd = 0;
      XMoveResizeWindow(d, ic[8], 320, 2, 40, 20);
      XUnmapWindow(d, w[8]);
    }
    break;
  }
}
/*--------------------------------------------------------------------*/

