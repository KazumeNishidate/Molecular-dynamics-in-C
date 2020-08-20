#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  initialize the system
*****/
void   init(void)
{
  get_control_param();  /* calculational control parameters: control.c */
  unit_converter();     /* parameter unit converter         */
  init_mem();           /* dynamic memory allocation        */
  set_potential();      /* potential setting: control.c     */
  identify_ion();       /* identify the ions: control.c     */
  set_loc();            /* initial locations: control.c     */
  set_vel();            /* initial velocities               */
  moment_correction();  /* velocity distribution correction */
}

/*****
  calculational unit conversion
*****/
void unit_converter(void)
/*****
*  purpose: unit conversion
*  all valiables are specified in the IEMD calculational unit.
*  see also [manual/units.IEMD]
*****/
{
  double press_unit_conversion;

  sys.Lx = sys.Ax * (double)sys.nx;   /* basic MD box size */
  sys.Ly = sys.Ay * (double)sys.ny;
  sys.Lz = sys.Az * (double)sys.nz;

  if(sys.a1 > 0) { /* if there is an ionic interaction [EWALD part] */
    sys.p2a2 = PIPI/( sys.a1 * sys.a1 ); 
    sys.a3   = sys.a1/SQRTPI;
  }

  sys.radius2 = sys.radius*sys.radius;
  sys.hm_sqrt = (short)sqrt((float)sys.hm); 

  sys.N =                 /* total number of ions           */
    ctl.natoms_in_unit_cell*sys.nx*sys.ny*sys.nz; 

  /* sys.perMol:  [energy] -> [J][mol^{-1}]  unit conversion    */
  sys.perMol = AVOGADRO/((double)(sys.N/ctl.natoms_in_mol_unit)); 

  sys.dt = ctl.delta_time_fs/SECD_2_FS;
  sys.dt2 = sys.dt*sys.dt;

  /* [GPa] * press_unit_conversion = [Pa']                             */
  press_unit_conversion = 1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );

  ctl.press_X = ctl.set_press_GPa_X*press_unit_conversion;
  ctl.press_Y = ctl.set_press_GPa_Y*press_unit_conversion;
  ctl.press_Z = ctl.set_press_GPa_Z*press_unit_conversion;  

  /* potential and force table for <Ewald 1st term> + <repulsion>    */
  /* 0.5 [A] upto sys.radius with 0.001[A] division.                  */
  sys.table_division = (int)((sys.radius - 0.5)*1000.0)+1;

  sys.kB = 1.380658e-5;   /* Boltzmann constant [IEMD unit] */
  sys.kk = 2.307079556;   /* e^2/(4 Pi epsilon) [IEMD unit] */

  ctl.kinds_of_interaction =   
    (ctl.kinds_of_ions*(ctl.kinds_of_ions+1))/2;

  /* [energy] -> [temperature] = [K] */
  /* T = (2/(3NKb))Sum(E_k), where Sum(E_k)=(1/2)Sum(mv^2) */
  sys.e2t = 2.0 / (3.0 * ((double)(sys.N) * sys.kB));

  /* [Pa'] -> [GPa] */
  sys.pp2gpa = (1.6605402*10000.0)/((SECD_2_FS)*(SECD_2_FS));

}

/*****
  function: void init_mem(void)
  purpose: dynamic memory allocation
*****/
void   init_mem(void)  
{
  /* to identify the atoms */
  sys.ion = (short *)calloc(sys.N+1,sizeof(short));

  /* mass of the atom */
  ion.m = (double *)calloc(ctl.kinds_of_ions, sizeof(double)); 
  /* charge of the atom */
  ion.Z = (double *)calloc(ctl.kinds_of_ions, sizeof(double)); 

  /* particle positions */	
  sys.rx = (double *)calloc(sys.N, sizeof(double)); 
  sys.ry = (double *)calloc(sys.N, sizeof(double)); 
  sys.rz = (double *)calloc(sys.N, sizeof(double)); 

  /* old forces */
  sys.fx0 = (double *)calloc(sys.N, sizeof(double)); 
  sys.fy0 = (double *)calloc(sys.N, sizeof(double)); 
  sys.fz0 = (double *)calloc(sys.N, sizeof(double)); 

  /* current forces */
  sys.fx = (double *)calloc(sys.N, sizeof(double)); 
  sys.fy = (double *)calloc(sys.N, sizeof(double)); 
  sys.fz = (double *)calloc(sys.N, sizeof(double)); 

  /* current velocities */
  sys.vx = (double *)calloc(sys.N, sizeof(double)); 
  sys.vy = (double *)calloc(sys.N, sizeof(double)); 
  sys.vz = (double *)calloc(sys.N, sizeof(double)); 

  /* for Ewald 2nd term calculation */	
  sys.rcsi = (double *)calloc(sys.N+1, sizeof(double)); 
  sys.rsni = (double *)calloc(sys.N+1, sizeof(double)); 


  /* potential and force table for <Ewald 1st term> + <repulsion>    */
  /* 0.5 [A] upto sys.radius with 0.001[A] division.                  */
  /*                                                                  */
  /* table_division = total number of divisions (column number)       */
  /* ddr = (int)[(0.5 [A] to sys.radius (cut-off radius)) - 0.5 [A]]  */

  /* look up table for potential */
  sys.pe1r1 = (double
	  *)calloc((ctl.kinds_of_interaction)*(sys.table_division+1),
		   sizeof(double));

  /* look up table for force */
  sys.fe1r1 = (double
	  *)calloc((ctl.kinds_of_interaction)*(sys.table_division+1),
		   sizeof(double));

  sys.lookup = (int
		*)calloc(ctl.kinds_of_ions*ctl.kinds_of_ions,
			 sizeof(int));

  /* for MSD calculation */
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
  set up initial velocities correspoding to the specified temperature
  value for each particles using the random real number sequence
  generator with normal (Gaussian) distribution.  see also "nrand()"
  in "program/ext.c".
*****/
void   set_vel(void)
{
  double bunsan = 2.0, sqrt_mass;
  int i;

  if(ctl.temp <= 0.0) return;

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

/***
 calc_alpha()
   purpose: estimation of alpha value for EWALD calculation
   see also the following function "reciprocal_space3a()" and
   "manual/manual.IEMD".
***/
void calc_alpha(void)
{
  double alpha=0.01;
  int i;

  for(i=0;i<100;i++)
    {
      alpha += 0.01;
      sys.a1   = alpha;
      sys.p2a2 = PIPI/( alpha * alpha ); 
      sys.a3   = alpha/SQRTPI;
      sys.pot = 0.0;
      
      mk_table();  

      real_space();            /* real space sum */
      reciprocal_space();      /* EWALD 2nd term */
      reciprocal_space3a();    /* EWALD 3rd term */

      printf("%f  %f\n",alpha,sys.pot);
    }
  exit(0);  /* EXIT the IEMD after the alpha estimation */
}

/***
 reciprocal_space3a()
   purpose: EWALD 3rd term calculation with changing the alpha value
   see also the above function "calc_alpha()"          
***/
void	reciprocal_space3a(void)
{
  static double	Z2, ewald_pot3 = 0.0;
  int i;

  Z2 = 0.0;
  for(i=0; i<sys.N; i++) {
    Z2 += ion_z(i) * ion_z(i);
  }
  ewald_pot3 = -Z2*sys.a3;  /* the value of sys.a3 changes with sys.a(alpha) */
  sys.pot += ewald_pot3;    /* potential contribution */
}





