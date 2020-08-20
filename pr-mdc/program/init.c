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
  init_mem1();          /* dynamic memory allocation 1                 */
  get_control_param();  /* calculational control parameters: control.c */
  unit_converter();     /* parameter unit converter                    */
  init_mem2();          /* dynamic memory allocation 2                 */

  set_potential();      /* potential setting: control.c                */
  identify_ion();       /* identify the ions: control.c                */
  set_loc();            /* initial locations: control.c                */

  evaluate_PR_param();  /* evaluate the cell mass and cell dt          */
  calc_cell();          /* evaluate the cell matrix : Hm -> Hmr        */

  set_vel();            /* initial velocities                          */
  moment_correction();  /* velocity distribution correction            */
}

/*****
  function: void init_mem1(void)
  purpose: dynamic memory allocation
*****/
void   init_mem1(void)
{
  pt.pres_set = dmat2d(3,3);     /* target pressure tensor              */
                                 /* 3x3 matrix : [0][0] to [2][2] array */
  pt.pres_tensor = dmat2d(3,3);  /* pressure tensor : Pab {a,b =x,y,z}  */
  pt.virial = dmat2d(3,3);       /* virial                              */
  pt.kin_tensor = dmat2d(3,3);   /* kinetic energy tensor               */

  cell.hu  = dmat2d(3,3);   /* unit cell matrix Hu (3x3) = {a,b,c} [A]  */
  cell.hur = dmat2d(3,3);   /* unit cell matrix Hu^(-1)                 */

  cell.hm      = dmat2d(3,3);   /* MD super cell matrix Hm              */
  cell.hmr     = dmat2d(3,3);   /* MD super cell matrix Hm^(-1)         */

  cell.hmt  = dmat2d(3,3);  /* transverse of Hm                         */
  cell.hmtr = dmat2d(3,3);  /* inverse of the transverse of Hm          */

  cell.hmv  = dmat2d(3,3);  /* transverse of (d/dt) Hm                  */
  cell.hmvt = dmat2d(3,3);  /* transverse of (d/dt) Hm                  */

  cell.hmf  = dmat2d(3,3);  /* cell force : [Hm^(T)]^(-1)*[(d/dt)Hm]^(T)*/

  cell.sigma = dmat2d(3,3); /* {sigma} matrix = { bxc, cxa, axb }       */
  cell.g     = dmat2d(3,3); /* {g} matrix = {sigma}/cell.vol            */

  /* [Parrinello Rahman method] */ 
  pr.v  = dmat2d(3,3);      /* cell velocity  */
  pr.f  = dmat2d(3,3);      /* cell force     */
  pr.f0 = dmat2d(3,3);      /* old cell force */

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
  int i, j;

  cell.hm[0][0] = cell.hu[0][0] * (double)cell.nx;   /* cubic MD super cell */
  cell.hm[1][1] = cell.hu[1][1] * (double)cell.ny;
  cell.hm[2][2] = cell.hu[2][2] * (double)cell.nz;

  sys.radius2 = sys.radius*sys.radius;
  sys.hm_sqrt = (int)sqrt((float)sys.hm); 

  sys.N =                 /* total number of ions           */
    ctl.natoms_in_unit_cell*cell.nx*cell.ny*cell.nz; 

  /* sys.perMol:  [energy] -> [J][mol^{-1}]  unit conversion    */
  sys.perMol = AVOGADRO/((double)(sys.N/ctl.natoms_in_mol_unit)); 

  sys.dt = ctl.delta_time_fs/SECD_2_FS;
  sys.dt2 = sys.dt*sys.dt;

  /* [GPa] * press_unit_conversion = [Pa']                            */
  press_unit_conversion = 1.0/(  (1.6605402*10000.0)/
                      ((SECD_2_FS)*(SECD_2_FS))  );

  for(i=0;i<3;i++) { /* setting pressure unit conversion */
    for(j=0;j<3;j++) {    
     pt.pres_set[i][j] *= press_unit_conversion;
    }
  }

  /* potential and force table for <repulsion>      */
  /* 0.5 [A] upto sys.radius with 0.001[A] division. */
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
  function: void init_mem2(void)
  purpose: dynamic memory allocation
*****/
void   init_mem2(void)  
{
  /* to identify the atoms */
  sys.ion = (int *)calloc(sys.N+1,sizeof(int));

  /* mass of the atom */
  ion.m = (double *)calloc(ctl.kinds_of_ions, sizeof(double)); 
  /* charge of the atom */
  ion.Z = (double *)calloc(ctl.kinds_of_ions, sizeof(double)); 

  /* particle positions */	
  cell.rx = (double *)calloc(sys.N, sizeof(double)); 
  cell.ry = (double *)calloc(sys.N, sizeof(double)); 
  cell.rz = (double *)calloc(sys.N, sizeof(double)); 

  cell.sx = (double *)calloc(sys.N, sizeof(double)); 
  cell.sy = (double *)calloc(sys.N, sizeof(double)); 
  cell.sz = (double *)calloc(sys.N, sizeof(double)); 

  /* current velocities */
  cell.vx = (double *)calloc(sys.N, sizeof(double)); 
  cell.vy = (double *)calloc(sys.N, sizeof(double)); 
  cell.vz = (double *)calloc(sys.N, sizeof(double)); 

  /* old forces */
  cell.fx0 = (double *)calloc(sys.N, sizeof(double)); 
  cell.fy0 = (double *)calloc(sys.N, sizeof(double)); 
  cell.fz0 = (double *)calloc(sys.N, sizeof(double)); 

  /* current forces */
  cell.fx = (double *)calloc(sys.N, sizeof(double)); 
  cell.fy = (double *)calloc(sys.N, sizeof(double)); 
  cell.fz = (double *)calloc(sys.N, sizeof(double)); 

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
  double bunsan = 5.0, sqrt_mass;
  int i;

  if(pt.temp_set_initial <= 0.0) return;

  for(i=0; i<sys.N; i++) {
    sqrt_mass = sqrt(ion_m(i));
    cell.vx[i] = nrand(pt.temp_set_initial, bunsan)/sqrt_mass;
    cell.vy[i] = nrand(pt.temp_set_initial, bunsan)/sqrt_mass;
    cell.vz[i] = nrand(pt.temp_set_initial, bunsan)/sqrt_mass;
  }
}

/*****
  clear old force
*****/
void   clear_foc(void)
{
  int i;
  for(i=0; i<sys.N; i++) {
    cell.fx0[i] = cell.fx[i];
    cell.fy0[i] = cell.fy[i]; 
    cell.fz0[i] = cell.fz[i]; 
    cell.fx[i] = 0.0;
    cell.fy[i] = 0.0;
    cell.fz[i] = 0.0;
  }
}

/*--------------------------------------------------------*/
/* memory allocation of 2-dimensional double format array */
/* set up the pointer for the pointer of array            */
/*--------------------------------------------------------*/
double **dmat2d(int n_row, int n_column)
{
  int i;
  double **mat;

  mat = (double **)calloc(n_row, sizeof(double*));
  if(!mat) { printf("calloc error in dmat2d\n"); exit(0); }
  
  mat[0] = (double *)calloc(n_row*n_column, sizeof(double));
  if(!mat[0]) { printf("calloc error in dmat2d\n"); exit(0); }

  for(i=1;i<n_row;i++) mat[i]=mat[i-1]+n_column;

  return mat;
}
