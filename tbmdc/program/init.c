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
  double *dmem, *ion_m_mem;
  short  *smem;

  dmem = (double *)calloc(12*sys.N+2, sizeof(double)); 
  smem = (short *)calloc(ctl.kinds_of_ions*(sys.N+1),sizeof(short));

  ion_m_mem = (double *)calloc(ctl.kinds_of_ions, sizeof(double)); 
  ion.m = ion_m_mem;

/* particle positions */	
  sys.rx = dmem;
  sys.ry = dmem + (sys.N);
  sys.rz = dmem + (sys.N) * 2;

/* old fources */
  sys.fx0 = dmem + (sys.N) * 3;
  sys.fy0 = dmem + (sys.N) * 4;
  sys.fz0 = dmem + (sys.N) * 5;

/* current fources */
  sys.fx = dmem + (sys.N) * 6;
  sys.fy = dmem + (sys.N) * 7;
  sys.fz = dmem + (sys.N) * 8;

/* current velocities */
  sys.vx = dmem + (sys.N) * 9;
  sys.vy = dmem + (sys.N) * 10;
  sys.vz = dmem + (sys.N) * 11;

  sys.ion = smem; 

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

  /* ----- taka (Monday, June 5, 2000) ----- */
  press.dx = (double *)calloc(sys.N*sys.N+1, sizeof(double));
  press.dy = (double *)calloc(sys.N*sys.N+1, sizeof(double));
  press.dz = (double *)calloc(sys.N*sys.N+1, sizeof(double));

  press.matX = (double *)calloc(sys.N*sys.N*16, sizeof( double ));
  press.matY = (double *)calloc(sys.N*sys.N*16, sizeof( double ));
  press.matZ = (double *)calloc(sys.N*sys.N*16, sizeof( double ));
  /*------------------------------------------*/
}

/*****
  set up initial velocities correcpoding to a specified temperature
  value for each particles using the random real number sequence
  generator of normal (Gaussian) distribution.  see also "nrand()"
  in "program/ext.c".
*****/
void   set_vel(void)
{
  double bunsan = 5.0, sqrt_mass;
  int i;

  for(i=0; i<sys.N; i++) {
    sqrt_mass = sqrt(ion_m(i));
    sys.vx[i] = nrand(ctl.temp, bunsan)/sqrt_mass;
    sys.vy[i] = nrand(ctl.temp, bunsan)/sqrt_mass;
    sys.vz[i] = nrand(ctl.temp, bunsan)/sqrt_mass;
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
  get_control_param(); /* calculational control parameters: control.c */
  init_mem();                     /* dynamic memory allocation        */
  set_potential();                /* potential setting: control.c     */
  identify_ion();                 /* identify the ions: control.c     */
  if(ctl.temp > 0.0) set_vel();   /*  initial velocities              */
  set_roc();                     /* initial locations: control.c     */
  /*  moment_correction();  */          /* velocity distribution correction */


}
