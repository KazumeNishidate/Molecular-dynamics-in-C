#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*********************************************************************/
/*  Carbon (diamond structure) with Abell Tersoff potential          */
/*********************************************************************/
/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*****/
void get_control_param(void)
{
  ctl.calc_max = 10000;       /* maximum MD time step                       */
  ctl.delta_time_fs = 0.38317;  /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 0.0;             /* Temperature setting [K]                    */
  ctl.t_control_step = 5;     /* scale at every ctl.t_control_step steps    */
  ctl.set_press_GPa_X = 0.01; /* Pressure setting  [GPa]                    */
  ctl.set_press_GPa_Y = 0.01; 
  ctl.set_press_GPa_Z = 0.01; 
  ctl.p_control_step  =  20;
  sys.pres_scalling_factor = 0.5;

  sys.Ax = 30.0;   /* initial lattice constant [C with c60 structure] */
  sys.Ay = 30.0;
  sys.Az = 30.0;
  sys.nx = 1;	   /* number of unit cells in x-direction                 */
  sys.ny = 1;      /* note: currently only nx=ny=nz=1 of C60 is supported */
  sys.nz = 1;

  ctl.natoms_in_unit_cell = 60;  /* number of atoms in unit cell          */
  ctl.natoms_in_mol_unit  = 60;  /* number of atoms in primitive Mol unit */
  ctl.kinds_of_ions       = 1;     /* C */

  /* alpha setting for EWALD calculation */
  sys.a1 = 0.0;     /* to SKIP the EWALD calculation */
  sys.hm = 20;      /* has no meaning for AT-C */
  sys.radius = 3.0; /* [A] AT potential set */  
}

/*****
  identify_ion:
  every atom is "C" in C crystal with diamond structure 
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
  set up C60 
*****/
void	set_loc(void)
{
  int i;
  /****
   *  C_60 structure
   *  ctl.natoms_in_unit_cell = C60
    ****/

  /* note: currently only nx=ny=nz=1 of C60 is supported */

  sys.rx[0]  = -2.59984; sys.ry[0]  = -2.33856; sys.rz[0]  = -0.58862; 
  sys.rx[1]  = -2.59984; sys.ry[1]  = -1.57232; sys.rz[1]  = -1.82843;
  sys.rx[2]  = -1.42071; sys.ry[2]  = -1.95544; sys.rz[2]  = -2.59468; 
  sys.rx[3]  = -0.69197; sys.ry[3]  = -2.95847; sys.rz[3]  = -1.82843; 
  sys.rx[4]  = -1.42071; sys.ry[4]  = -3.19525; sys.rz[4]  = -0.58862; 
  sys.rx[5]  =  2.59984; sys.ry[5]  = -2.33856; sys.rz[5]  = -0.58862; 
  sys.rx[6]  =  1.42071; sys.ry[6]  = -3.19525; sys.rz[6]  = -0.58862;
  sys.rx[7]  =  0.69197; sys.ry[7]  = -2.95847; sys.rz[7]  = -1.82843;
  sys.rx[8]  =  1.42071; sys.ry[8]  = -1.95544; sys.rz[8]  = -2.59468;
  sys.rx[9]  =  2.59984; sys.ry[9]  = -1.57232; sys.rz[9]  = -1.82843;
  sys.rx[10] = -3.47789; sys.ry[10] =  0.36379; sys.rz[10] = -0.58862;
  sys.rx[11] = -3.02750; sys.ry[11] =  1.74994; sys.rz[11] = -0.58862;
  sys.rx[12] = -2.29876; sys.ry[12] =  1.98672; sys.rz[12] = -1.82843;
  sys.rx[13] = -2.29876; sys.ry[13] =  0.74691; sys.rz[13] = -2.59468;
  sys.rx[14] = -3.02750; sys.ry[14] = -0.25612; sys.rz[14] = -1.82843;
  sys.rx[15] =  3.47789; sys.ry[15] =  0.36379; sys.rz[15] = -0.58862;
  sys.rx[16] =  3.02750; sys.ry[16] = -0.25612; sys.rz[16] = -1.82843;
  sys.rx[17] =  2.29876; sys.ry[17] =  0.74691; sys.rz[17] = -2.59468;
  sys.rx[18] =  2.29876; sys.ry[18] =  1.98672; sys.rz[18] = -1.82843;
  sys.rx[19] =  3.02750; sys.ry[19] =  1.74994; sys.rz[19] = -0.58862;
  sys.rx[20] =  3.47789; sys.ry[20] = -0.36379; sys.rz[20] =  0.58862;
  sys.rx[21] =  3.02750; sys.ry[21] =  0.25612; sys.rz[21] =  1.82843;
  sys.rx[22] =  2.29876; sys.ry[22] = -0.74691; sys.rz[22] =  2.59468;
  sys.rx[23] =  2.29876; sys.ry[23] = -1.98672; sys.rz[23] =  1.82843;
  sys.rx[24] =  3.02750; sys.ry[24] = -1.74994; sys.rz[24] =  0.58862;
  sys.rx[25] = -3.02750; sys.ry[25] = -1.74994; sys.rz[25] =  0.58862;
  sys.rx[26] = -2.29876; sys.ry[26] = -1.98672; sys.rz[26] =  1.82843;
  sys.rx[27] = -2.29876; sys.ry[27] = -0.74691; sys.rz[27] =  2.59468;
  sys.rx[28] = -3.02750; sys.ry[28] =  0.25612; sys.rz[28] =  1.82843;
  sys.rx[29] = -3.47789; sys.ry[29] = -0.36379; sys.rz[29] =  0.58862;
  sys.rx[30] = -0.72874; sys.ry[30] = -1.00303; sys.rz[30] = -3.32225;
  sys.rx[31] = -1.17913; sys.ry[31] =  0.38312; sys.rz[31] = -3.32225;
  sys.rx[32] = -0.00000; sys.ry[32] =  1.23981; sys.rz[32] = -3.32225;
  sys.rx[33] =  1.17913; sys.ry[33] =  0.38312; sys.rz[33] = -3.32225;
  sys.rx[34] =  0.72874; sys.ry[34] = -1.00303; sys.rz[34] = -3.32225;
  sys.rx[35] = -0.72874; sys.ry[35] = -3.42008; sys.rz[35] =  0.58862;
  sys.rx[36] =  0.72874; sys.ry[36] = -3.42008; sys.rz[36] =  0.58862;
  sys.rx[37] =  1.17913; sys.ry[37] = -2.80018; sys.rz[37] =  1.82843;
  sys.rx[38] = -0.00000; sys.ry[38] = -2.41705; sys.rz[38] =  2.59468;
  sys.rx[39] = -1.17913; sys.ry[39] = -2.80018; sys.rz[39] =  1.82843;
  sys.rx[40] =  1.17913; sys.ry[40] = -0.38312; sys.rz[40] =  3.32225;
  sys.rx[41] =  0.72874; sys.ry[41] =  1.00303; sys.rz[41] =  3.32225;
  sys.rx[42] = -0.72874; sys.ry[42] =  1.00303; sys.rz[42] =  3.32225;
  sys.rx[43] = -1.17913; sys.ry[43] = -0.38312; sys.rz[43] =  3.32225;
  sys.rx[44] =  0.00000; sys.ry[44] = -1.23981; sys.rz[44] =  3.32225;
  sys.rx[45] = -0.72874; sys.ry[45] =  3.42008; sys.rz[45] = -0.58862;
  sys.rx[46] =  0.72874; sys.ry[46] =  3.42008; sys.rz[46] = -0.58862;
  sys.rx[47] =  1.17913; sys.ry[47] =  2.80018; sys.rz[47] = -1.82843;
  sys.rx[48] = -0.00000; sys.ry[48] =  2.41705; sys.rz[48] = -2.59468;
  sys.rx[49] = -1.17913; sys.ry[49] =  2.80018; sys.rz[49] = -1.82843;
  sys.rx[50] = -1.42071; sys.ry[50] =  1.95544; sys.rz[50] =  2.59468;
  sys.rx[51] = -0.69197; sys.ry[51] =  2.95847; sys.rz[51] =  1.82843;
  sys.rx[52] = -1.42071; sys.ry[52] =  3.19525; sys.rz[52] =  0.58862;
  sys.rx[53] = -2.59984; sys.ry[53] =  2.33856; sys.rz[53] =  0.58862;
  sys.rx[54] = -2.59984; sys.ry[54] =  1.57232; sys.rz[54] =  1.82843;
  sys.rx[55] =  2.59984; sys.ry[55] =  2.33856; sys.rz[55] =  0.58862;
  sys.rx[56] =  1.42071; sys.ry[56] =  3.19525; sys.rz[56] =  0.58862;
  sys.rx[57] =  0.69197; sys.ry[57] =  2.95847; sys.rz[57] =  1.82843;
  sys.rx[58] =  1.42071; sys.ry[58] =  1.95544; sys.rz[58] =  2.59468;
  sys.rx[59] =  2.59984; sys.ry[59] =  1.57232; sys.rz[59] =  1.82843;
  
  for(i=0; i<ctl.natoms_in_unit_cell; i++){ 
    sys.rx[i] += sys.Ax/2.0; /* set the origin to the center of cell */
    sys.ry[i] += sys.Ay/2.0;
    sys.rz[i] += sys.Az/2.0;
  }

}

void set_potential(void)
{
  double m_C;

  m_C = 12.011 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */

  ion.m[0] = m_C;

  /* ------------ A.T Potential Set --------------*/
  at.A = 1.3936e3 * 1.60217733e-1;  /* [eV -> en] */
  at.B = 3.4670e2 * 1.60217733e-1;  /* [eV -> en] */
  at.lam = 3.4879;
  at.mu = 2.2119;
  at.beta = 1.5724e-7;
  at.n = 7.2751e-1;
  at.c = 3.8049e4;
  at.d = 4.384;
  at.h = -0.57058;
  at.R = 1.8;
  at.S = 2.1;               /* = sys.radius = cut-off radius */
  at.x = 1.0;
  /* ---------------------------------------------*/

}

void mk_table(void)
{
  /*  Carbon (diamond structure) with Abell Tersoff potential      */
  /*  this function do nothing.  see "real.c".                     */
}






