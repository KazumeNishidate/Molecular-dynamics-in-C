#include  <X11/Xlib.h>
#include  <X11/Xutil.h>
#include  <stdio.h>
#include  "../prototypes.h"

unsigned long MyColor(xdisplay, color)
     Display *xdisplay;
     char *color;
{
  Colormap cmap;
  XColor c0, c1;
  
  cmap = DefaultColormap (xdisplay, 0);
  XAllocNamedColor (xdisplay, cmap, color, &c1, &c0); 
  return(c1.pixel);
}
