#include  <stdio.h>
#include  <sys/types.h>
#include  <time.h>
#include  "potential.h"

/*------------------ global definitions ------------------------------------*/
#ifndef  PI
#define  PI         3.1415926535897932385
#endif
#ifndef  PI2
#define  PI2        6.283185307179587       /* = 2*Pi                      */
#endif
#define  PIPI       9.869604401089358619    /* Pi*Pi                       */
#define  SQRTPI     1.77245385090551602731  /* Sqrt(Pi)                    */
#define  SECD_2_FS  4.074972637944947       /* second' (IEMD unit) -> fsec */

#define  AVOGADRO   6.0221367e5             /* Avogadro constant in IEMD unit  */

/* input and output files */
FILE *fpout, *fpout_positions, *fpout_velocities, *fpmac;

/*------------------ global variables --------------------------------------*/
typedef	struct {

  double delta_time_fs;     /* delta_t [fsec]                          */
  int calc_max;             /* maximum MD step number                  */
  double temp;              /* temperature setting [K]                 */
  int t_control_step;  
  int p_control_step;

  double set_press_GPa_X;   /* pressure setting [GPa]                  */
  double set_press_GPa_Y;
  double set_press_GPa_Z;

  double press_X;           /* pressure setting [TPa] (IEMD unit)      */
  double press_Y;
  double press_Z;

  int kinds_of_ions;        /* total number of kinds of ions           */

  int natoms_in_mol_unit;   /* number of atoms in primitive unit cell  */
                            /* [to calculate the (Mol) unit]           */

  int natoms_in_unit_cell;  /* number of atoms in unit cell            */
               /* sys.N = ctl.natoms_in_unit_cell*sys.nx*sys.ny*sys.nz */

  /* In a system of par-wise potential,                                */
  /*   kinds_of_interaction = Permutation_{MAX_IONS, MAX_IONS}         */
  int kinds_of_interaction; /* total number of kinds of interaction    */

} calc_control;

typedef	struct {
  double  *rx, *ry, *rz;    /* positions                           */
  double  *vx, *vy, *vz;    /* velocities                          */
  double  *fx0, *fy0, *fz0; /* old force                           */
  double  *fx, *fy, *fz;    /* current foce                        */
  
  short	*ion;	 /* kind of ions                                   */
  short	N;	 /* total number of ions in the system             */
  short	nx;	 /* number of unit cells in x-direction            */
  short	ny;	 /*                         y-direction            */
  short	nz;	 /*                         z-direction            */

  double  Ax;	 /* initial lattice constant of x-direction        */
  double  Ay;	 /*                             y-direction        */
  double  Az;	 /*                             z-direction        */
  double  Lx;	 /* basic MD cell size = sys.Ax * (double)sys.nx   */
  double  Ly;	 /*                      sys.Ay * (double)sys.ny   */
  double  Lz;	 /*                      sys.Az * (double)sys.nz   */

  double  kB;    /* Boltzman constant                              */
  double  kk;    /* kk = e^2/(4Pi*Epsiron) [m/F]  <MKS-unit>       */
                 /*    = 1 [cm*erg/(stat.C^2)]  <cgs-unit>         */

  double  perMol;  /* [energy] -> [J][mol^{-1}]                    */

  double  e2t;     /* [energy] -> [temperature] = [K] converter    */
  double  e2tX;
  double  e2tY;
  double  e2tZ;
  double  pp2gpa;  /* pressure unit conversion : [TPa] -> [GPa]    */

  int     step;    /* store the current MD time step number        */
  double  dt;	   /* delta time                                   */
  double  dt2;	   /* delta time^2                                 */

  double  a1;      /* alpha parameter for EWALD calculation        */
  double  p2a2;    /* = Pi^2/alpha^2 for EWALD calculation         */
  double  a3;      /* alpha/sqrt(pi)                               */  

  short  hm;       /* = (maximum size of reciprocal vector)^2      */
  short  hm_sqrt;  /* = sqrt(maximum size of reciprocal vector)    */

  double  radius;  /* potential cutoff radious for real-space      */
  double  radius2; /* radius^2                                     */

  /* look-up table of the reciprocal space calculation in EWALD    */
  double  *rcsi, *rsni;     

  /* total number of divisions (column number);                     */
  /* potential and force table for <Ewald 1st term> + <repulsion>  */
  /* of r = (0.5 [A] upto sys.radius) region with 0.001[A] division.*/
  /* see also init_mem(), mk_table(), real_space()                  */
  int  table_division; 

  /* table for potential calculation which will be created by      */
  /* mk_table1().  <Ewald-first-term> + <repulsion>                */
  /* see also init.c real.c, md.h                                  */
  double  *pe1r1; 

  /* table for force calculation which will be created by          */
  /* mk_table1(). <Ewald-first-term> + <repulsion>                 */
  /* see also init.c real.c, md.h                                  */
  double  *fe1r1;  

  /* lookup-table index */
  int *lookup;

  double  virX;  /* virial */
  double  virY;
  double  virZ;

  double  presX; /* pressure tensor P_xx */
  double  presY; /*                 P_yy */
  double  presZ; /*                 P_zz */
  double  pres;  /* current pressure     */

  /* cell size scalling factor referring the stress tensors Pxx, Pyy,
     and Pzz */
  double  pres_scalling_factor;   

  double  pot;
  double  kin;
  double  kinX;
  double  kinY;
  double  kinZ;

} system_property;

typedef	struct {
  double	*m;	/* mass                                    */
  double	*Z;	/* effective charge                        */
} ion_property;

typedef	struct {
  short	*pos;          /* positions                                */
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
  double Lx2, Ly2, Lz2;
} msd_calculation_variables;

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
  ion_property ion;
  system_property sys;
  calc_control ctl;
  network_resolver net;
  msd_calculation_variables msd;
  Gear7F gear;

/*------------------- macro -----------------------------------------------*/
#define	ion_m(i)	(ion.m[sys.ion[i]])
#define	ion_z(i)	(ion.Z[sys.ion[i]])
/*-------------------------------------------------------------------------*/
