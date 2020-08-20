#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  function: void init_mem(void)
  purpose: dynamic memory allocation
*****/
void   init_mem(void)  
{

  ion.m = (double *)calloc(ctl.kinds_of_ions, sizeof(double)); 

/* particle positions */	
  sys.rx = (double *)calloc(sys.N, sizeof(double)); 
  sys.ry = (double *)calloc(sys.N, sizeof(double)); 
  sys.rz = (double *)calloc(sys.N, sizeof(double)); 

/* old fources */
  sys.fx0 = (double *)calloc(sys.N, sizeof(double)); 
  sys.fy0 = (double *)calloc(sys.N, sizeof(double)); 
  sys.fz0 = (double *)calloc(sys.N, sizeof(double)); 

/* current fources */
  sys.fx = (double *)calloc(sys.N, sizeof(double)); 
  sys.fy = (double *)calloc(sys.N, sizeof(double)); 
  sys.fz = (double *)calloc(sys.N, sizeof(double)); 

/* current velocities */
  sys.vx = (double *)calloc(sys.N, sizeof(double)); 
  sys.vy = (double *)calloc(sys.N, sizeof(double)); 
  sys.vz = (double *)calloc(sys.N, sizeof(double)); 

  sys.ion = (short *)calloc(ctl.kinds_of_ions*(sys.N+1),sizeof(short));

  /* for MSD calculation */
  /* initial positions at 0 [step] (time=0) */
  msd.rx0 = (double *)calloc(sys.N+1, sizeof(double));     
  msd.ry0 = (double *)calloc(sys.N+1, sizeof(double));     
  msd.rz0 = (double *)calloc(sys.N+1, sizeof(double));     

  /* old positions at t-1 [step] (time=t-1) */
  msd.rx_old = (double *)calloc(sys.N+1, sizeof(double));     
  msd.ry_old = (double *)calloc(sys.N+1, sizeof(double));     
  msd.rz_old = (double *)calloc(sys.N+1, sizeof(double));     
  
  msd.dx = (double *)calloc(sys.N+1, sizeof(double));     
  msd.dy = (double *)calloc(sys.N+1, sizeof(double));     
  msd.dz = (double *)calloc(sys.N+1, sizeof(double));     
  msd.value = (double *)calloc(ctl.kinds_of_ions, sizeof(double));     
  msd.number_of_ion = (int *)calloc(ctl.kinds_of_ions, sizeof(int));     

}

/*****
  set up initial velocities correcpoding to a specified temperature
  value for each particles using the random real number sequence
  generator of normal (Gaussian) distribution.  see also "nrand()"
  in "program/ext.c".
*****/

void   set_vel(void)
{
  double nrand(double, double);
  double bunsan = 5.0;
  double mavrg = 0.0 ;
  int i;

  for(i=0; i<sys.N; i++) {
    mavrg += sqrt(ion_m(i));
  }
  mavrg /= sys.N;  /* = (Sum[sqrt(mass)]/sys.N) */

/* BUG */

  for(i=0; i<sys.N; i++) {
    sys.vx[i] = nrand(ctl.temp, bunsan)/mavrg;
    sys.vy[i] = nrand(ctl.temp, bunsan)/mavrg;
    sys.vz[i] = nrand(ctl.temp, bunsan)/mavrg;
  }

}

/*****
  clear old force
*****/
void   clear_foc(void)
{
  short i;
  for(i=0; i<sys.N; i++) {
    sys.fx0[i] = sys.fx[i];
    sys.fy0[i] = sys.fy[i]; 
    sys.fz0[i] = sys.fz[i]; 
    sys.fx[i] = 0.0;
    sys.fy[i] = 0.0;
    sys.fz[i] = 0.0;
  }
}

/*****
  initialize the system
*****/
void   init(void)
{
  get_control_param();  /* calculational control parameters: control.c */
  init_mem();           /* dynamic memory allocation        */
  set_potential();      /* potential setting: control.c     */
  identify_ion();       /* identify the ions: control.c     */
  /*  set_vel(); */           /* initial velocities               */
  set_roc();            /* initial locations: control.c     */
  moment_correction();  /* velocity distribution correction */

}
