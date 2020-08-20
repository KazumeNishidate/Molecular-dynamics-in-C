/*--- XD root window ---*/
#define  XD_MAXPTS 10000  /* maximum number of ions to display in XD */
#define  XD_NETS_IONS_COLOR 16
#define  XD_UNIT_IONS_COLOR 72
#define  XD_ROOT_WINDOWS 9
#define  XD_ROOT_ICONS 12
#define  XD_XNET_ICONS 12
#define  XD_XUNIT_ICONS 11
#define  XD_MAX_CHAR_NUM 100
#define  GRADATION 7       /* color gradation degree to use in XUNIT */

/*--- XHIST ---*/
#define EK_MAX 100          /* set Kinetic Energy scale Max [KJ/mol] */
#define EP_MAX -2500        /* set Kinetic Energy scale Max [KJ/mol] */
#define ET_MAX -2500        /* set Kinetic Energy scale Max [KJ/mol] */
#define TEMP_MAX 2000       /* set temperature scale Max [K] */
#define PRES_MAX 5          /* set pressure scale Max [GPa]  */

/*------------  X Window related declaration  ---------------------------*/
 Display *d;
 GC gc;
 XEvent e;
 XSetWindowAttributes attr;
 Window w[XD_ROOT_WINDOWS], ic[XD_ROOT_ICONS]; 
 Window ic_xnet[XD_XNET_ICONS], ic_xunit[XD_XUNIT_ICONS];
 char stringbuffer[XD_MAX_CHAR_NUM];
 unsigned long color[XD_NETS_IONS_COLOR];
 XColor   Col[XD_UNIT_IONS_COLOR];
 Colormap cmap;

/*----------------- other global variables ------------------------------*/
/*--- XNET ---*/
 double xnet_theta, xnet_phy;  /* xdinit.c xdicon.c xnet.c     */
 double xnet_zoom;             /* xdinit.c xdicon.c xnet.c       */
 int xnet_shift_x, xnet_shift_y;

/*--- XUNIT ---*/
 double xunit_theta, xunit_phy;  /* xdinit.c xdicon.c xunit.c     */
 double xunit_zoom;             /* xdinit.c xdicon.c xunit.c       */
 int xunit_shift_x, xunit_shift_y;

/*----------------- other global variables ------------------------------*/
/*  "enum" declaration to use in XNET                                      */
/*  Change the following order to select a color for each of ions in XNET  */
enum{Red, Blue, Green, Yellow, Orange, Cyan, BlueViolet, LimeGreen,
	Coral, Khaki, SpringGreen, LightGray, Magenta, SlateBlue, Black,
	White};
/*-----------------------------------------------------------------------*/

