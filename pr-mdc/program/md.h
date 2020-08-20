#include  <stdio.h>
#include  <sys/types.h>
#include  <time.h>
#include  "potential.h"

/*------------------- global definitions ------------------------------------*/
#ifndef  PI
#define  PI         3.1415926535897932385
#endif
#ifndef  PI2
#define  PI2        6.283185307179587      /* = 2*Pi                         */
#endif
#define  PIPI       9.869604401089358619   /* Pi*Pi                          */
#define  SQRTPI     1.77245385090551602731 /* Sqrt(Pi)                       */
#define  SECD_2_FS  4.074972637944947      /* second' (IEMD unit) -> fsec    */

#define  AVOGADRO   6.0221367e5            /* Avogadro constant in IEMD unit */

#define  RADIAN_TO_DEGREE(x) ((x)*180.0/PI)

/* input and output files */
FILE *fpout, *fpout_positions, *fpout_velocities, *fpout_cell;

/*------------------- global variables --------------------------------------*/
typedef	struct {

  int skip_ewald;           /* 0 => calc Ewald, 1 => skip Ewald        */
  int calc_max;             /* maximum MD step number                  */
  double delta_time_fs;     /* delta_t [fsec]                          */

  int kinds_of_ions;        /* total number of kinds of ions           */

  int natoms_in_mol_unit;   /* number of atoms in primitive unit cell  */
                            /* [to calculate the (Mol) unit]           */

  int natoms_in_unit_cell;  /* number of atoms in unit cell            */
               /* sys.N = ctl.natoms_in_unit_cell*sys.nx*sys.ny*sys.nz */

  /* In a system of par-wise potential,                                */
  /*   kinds_of_interaction = Permutation_{MAX_IONS, MAX_IONS}         */
  int kinds_of_interaction; /* total number of kinds of interaction    */

} calc_control;

typedef struct {
  int	nx;	 /* number of unit cells in x-direction             */
  int	ny;	 /*                         y-direction             */
  int	nz;	 /*                         z-direction             */

  double  vol;   /* MD cell volume                                  */

  double  *rx, *ry, *rz;    /* positions                            */
  double  *sx, *sy, *sz;    /* scalled positions                    */
  double  *vx, *vy, *vz;    /* velocities                           */
  double  *fx0, *fy0, *fz0; /* old force                            */
  double  *fx, *fy, *fz;    /* current foce                         */

  /*=== Cartesian coordinate{r} <-> scaled coordinate{s}          ===*/
  double **hu;    /* unit cell matrix Hu (3x3) = {a,b,c} [A]         */
                  /* initial setting = cubic unit cell               */
  double **hur;   /* inverse of the unit cell matrix [Hu^(-1)]       */

  double **hm;    /* MD super cell matrix Hm                         */
  double **hmr;   /* inverse of the MD super cell matrix [Hm^(-1)]   */

  double **hmt;   /* Hm^T : transverse of the MD super cell matrix   */
  double **hmtr;  /* (Hm^T)^(-1) : inverse of the Hm^T               */

  double **hmv;   /* (d/dt) Hm                                       */
  double **hmvt;  /* transverse of (d/dt) Hm                         */

  double **hmf;   /* cell force : [Hm^(T)]^(-1)*[(d/dt)Hm]^(T)       */

  double **g;     /* {g} matrix     : g ={bxc, cxa, axb}/omega       */
  double **sigma; /* {sigma} matrix : sigma ={bxc, cxa, axb}         */

  double Lx;
  double Ly;
  double Lz;

  double angle_ab;
  double angle_bc;
  double angle_ca;

} cell_parameters;

typedef	struct {
  int	*ion;	   /* kind of ions                                  */
  int	N;	   /* total number of ions in the system            */

  int     step;    /* store the current MD time step number         */
  double  dt;	   /* delta time                                    */
  double  dt2;	   /* delta time^2                                  */

  double  kB;      /* Boltzman constant                             */
  double  kk;      /* kk = e^2/(4Pi*Epsiron) [m/F]  <MKS-unit>      */
                   /*    = 1 [cm*erg/(stat.C^2)]  <cgs-unit>        */

  double  pot;     /* total potential                               */
  double  kin;     /* kinetic energy = Trace[pt.kin_tensor]         */

  double  perMol;  /* [energy] -> [J][mol^{-1}]                     */

  double  e2t;     /* [energy] -> [temperature] = [K] converter     */
  double  pp2gpa;  /* pressure unit conversion : [Pa'] -> [GPa]     */

  double  a1;      /* alpha parameter for EWALD calculation         */

  int  hm;         /* = (maximum size of reciprocal vector)^2       */
  int  hm_sqrt;    /* = sqrt(maximum size of reciprocal vector)     */

  double  radius;  /* cutoff radius for real-space table            */
  double  radius2; /* radius^2                                      */

  double  cutoff;  /* cutoff radius for                             */
  double  cutoff2; /*        short range repulsion + EWALD 1st term */

  /* look-up table of the reciprocal space calculation in EWALD     */
  double  *rcsi, *rsni;     

  /* total number of divisions (column number);                     */
  /* potential and force table for <Ewald 1st term> + <repulsion>  */
  /* of r = (0.5 [A] upto sys.radius) region with 0.001[A] division.*/
  /* see also init_mem(), mk_table(), real_space()                  */
  int  table_division; 

  /* table for potential calculation which will be created by       */
  /* mk_table1().  <Ewald-first-term> + <repulsion>                 */
  /* see also init.c real.c, md.h                                   */
  double  *pe1r1; 

  /* table for force calculation which will be created by           */
  /* mk_table1(). <Ewald-first-term> + <repulsion>                  */
  /* see also init.c real.c, md.h                                   */
  double  *fe1r1;  

  /* lookup-table index */
  int *lookup;

} system_property;

typedef struct {
  int t_control_step;  
  int p_control_step;

  double  pres_scalling_factor;  /* cell size scalling factor */

  double temp_set;          /* temperature setting [K]                 */
  double temp_set_initial;  /* temperature setting [K]                 */

  double pres_trace;        /* pressure = Trace[pressure tensor]       */
  double **pres_set;        /* pressure setting [GPa'] : 3x3 matrix    */
  double **pres_tensor;     /* pressure tensor  [GPa'] : 3x3 matrix    */
  double **virial;          /* virial */

  double  **kin_tensor;     /* kinetic energy tensor                   */

} pressure_temperature_set;

typedef	struct {
  double	*m;	/* mass              */
  double	*Z;	/* effective charge  */
} ion_property;

typedef	struct {
  int	  *pos;        /* positions                                */
  int     *x_pos;      /* scalled x-y positions to use in xroll.c  */
  int     *y_pos;
  int     counter;
  double  length;
} network_resolver;

typedef	struct {
  /* --------------------------------------------------------------- */
  /* Mean Squared Displacement calculation                           */
  /*                                                                 */
  /* MSD(t) = <|r(t)-r(0)|^2> => 6 D t   [Einstein's relation]       */
  /*                                                                 */
  /* diffusion constant D = <|r(t)-r(0)|^2>/(6t)                     */
  /* see also "calc_msd()" in "rv.c"                                 */
  /*                                                                 */
  /* note: this is a optional calculation. to calculate MSD, just    */
  /* comment out the corresponding function-call line from "newton()"*/
  /* in "main.c".                                                    */
  /* --------------------------------------------------------------- */
  int *number_of_ion;
  double *rx_old, *ry_old, *rz_old;
  double *dx, *dy, *dz;
  double *value;
} msd_calculation_variables;

typedef	struct {  /* pressure control by Parrinello Rahman method */
  double mass;    /* fictitious cell mass                   */
  double dt;      /* fictitious cell dt                     */
  double **v;     /* fictitious cell velocity  : 3x3 matrix */
  double **f;     /* fictitious cell force     : 3x3 matrix */
  double **f0;    /* fictitious old cell force : 3x3 matrix */
} Parrinello_Rahman_set;

typedef	struct {  /* to use in the Gear's algorithm */
  double  *anx;
  double  *any;
  double  *anz;

  double  *an1x;
  double  *an1y;
  double  *an1z;

  double  *an2x;
  double  *an2y;
  double  *an2z;

  double  *an3x;
  double  *an3y;
  double  *an3z;

  double  *an4x;
  double  *an4y;
  double  *an4z;
} Gear7F;

/*------------------- declaration for the structures ----------------------*/

  calc_control ctl;
  cell_parameters cell;

  system_property sys;
  pressure_temperature_set pt;
  ion_property ion;

  network_resolver net;
  msd_calculation_variables msd;

  Parrinello_Rahman_set pr;

  Gear7F gear;

/*------------------- macro -----------------------------------------------*/
#define	ion_m(i)	(ion.m[sys.ion[i]])
#define	ion_z(i)	(ion.Z[sys.ion[i]])
/*-------------------------------------------------------------------------*/
