#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*********************************************************************/
/*  Silicon (diamond structure) with Abell Tersoff potential         */
/*********************************************************************/
/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*****/
void get_control_param(void)
{
  ctl.calc_max = 10000;        /* maximum MD time step                       */
  ctl.delta_time_fs = 0.38317; /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 100.0;           /* Temperature setting [K]                    */
  ctl.t_control_step = 5;     /* scale at every ctl.t_control_step steps    */
  ctl.set_press_GPa_X = 0.01; /* Pressure setting  [GPa]                    */
  ctl.set_press_GPa_Y = 0.01; 
  ctl.set_press_GPa_Z = 0.01; 
  ctl.p_control_step  =  20;
  sys.pres_scalling_factor = 0.5;

  sys.Ax = 5.42;   /* initial lattice constant [Si with diamond structure] */
  sys.Ay = 5.42;
  sys.Az = 5.42;
  sys.nx = 3;	   /* number of unit cells in x-direction */
  sys.ny = 3;
  sys.nz = 3;

  ctl.natoms_in_unit_cell = 8;     /* number of atoms in unit cell          */
  ctl.natoms_in_mol_unit  = 2;     /* number of atoms in primitive Mol unit */
  ctl.kinds_of_ions = 1;   /* "Si" : see also "identify_ion()" */

  /* alpha setting for EWALD calculation */
  sys.a1   = 0.0;   /* to SKIP the EWALD calculation */
  sys.hm = 20;      /* has no meaning for AT-Si */
  sys.radius = 3.0; /* [A] AT potential set */
}

/*****
  identify_ion:
  every atom is "Si" in Si crystal with diamond structure 
*****/
void	identify_ion(void)
{
  short	x, y, z, i=0;

  for(z=0; z<sys.nz*2; z++)
    for(y=0; y<sys.ny*2; y++)
      for(x=0; x<sys.nx*2; x++)
	sys.ion[i++] = 0;
}

/*****
  void set_loc(void)
  set up a diamond structure for Si 
*****/
void	set_loc(void)
{
  short	x, y, z, i=0;
  short nx, ny, nz;
  double nnx=0.0, nny=0.0, nnz=0.0;
  double nx3=0.0, ny3=0.0, nz3=0.0;
  double shiftx, shifty, shiftz;

  nx = sys.nx * 2;
  ny = sys.ny * 2;
  nz = sys.nz * 2;

/**** 
  X display shift 
****/
  shiftx = (sys.Ax)/2.0;
  shifty = (sys.Ay)/2.0;
  shiftz = (sys.Az)/2.0;

/****
  fcc ion
****/

  for(z=0; z<nz; z++)
    for(y=0; y<ny; y++)
      for(x=0; x<nx; x++)
	{
	  if((x+y+z)%2==0){
	    nnx = (double)(x)/2.0;
	    nny = (double)(y)/2.0;
	    nnz = (double)(z)/2.0;
	    sys.rx[i] = sys.Ax * nnx +shiftx;
	    sys.ry[i] = sys.Ay * nny +shifty;
	    sys.rz[i] = sys.Az * nnz +shiftz;
	    i++;

	    nx3 = nnx + 0.25;
	    ny3 = nny + 0.25;
	    nz3 = nnz + 0.25;
	    sys.rx[i] = sys.Ax * nx3 +shiftx;
	    sys.ry[i] = sys.Ay * ny3 +shifty;
	    sys.rz[i] = sys.Az * nz3 +shiftz;
	    i++;
	  }
	}
}

void set_potential(void)
{
  double m_C;

  m_C = 28.09 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */

  ion.m[0] = m_C;

  /* ------------ A.T Potential Set --------------*/
  at.A = 1.8308e3 * 1.60217733e-1;  /* [eV -> en] */
  at.B = 4.7118e2 * 1.60217733e-1;  /* [eV -> en] */
  at.lam = 2.4799;
  at.mu = 1.7322;
  at.beta = 1.10e-6;          
  at.n = 7.8734e-1;
  at.c = 1.0039e5;
  at.d = 16.217;
  at.h = -0.59825;
  at.R = 2.7;            
  at.S = 3.0;               /* = sys.radius = cut-off radius */
  at.x = 1.0;

  /* ---------------------------------------------*/

}

void mk_table(void)
{
  /*  Silicon (diamond structure) with Abell Tersoff potential     */
  /*  this function do nothing.  see "real.c".                     */
}






