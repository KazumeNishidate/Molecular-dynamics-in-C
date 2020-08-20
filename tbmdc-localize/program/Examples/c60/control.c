#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*********************************************************************/
/* TBMD for Carbon                                                   */
/* [Ref.]                                                            */
/*        C. H. Xu, C. Z. Wang, C. T. Chan, and K. M. Ho,            */
/*        ``A transferable tight-binding potential for carbon'',     */
/*        J. Phys. Condens. Matter, Vol. 4, 6047 (1992).             */
/*********************************************************************/
/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*****/
void get_control_param(void)
{
  double set_press_GPa_X, set_press_GPa_Y, set_press_GPa_Z;
  double press_unit_conversion;
  double delta_time_fs;

  ctl.calc_max = 100000;  /* maximum MD time step                     */
  delta_time_fs = 0.35;   /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 300.0;       /* Temperature setting [K] */
  ctl.t_control_step = 5;
  set_press_GPa_X = 0.01;  /* Pressure setting  [GPa] */
  set_press_GPa_Y = 0.01; 
  set_press_GPa_Z = 0.01; 
  ctl.p_control_step = 20;
  sys.pres_scalling_factor = 0.5;

  /* check */
  sys.Ax = 30.0;   /* initial lattice constant */
  sys.Ay = 30.0;     /* initial lattice constant */
  sys.Az = 30.0;     /* initial lattice constant */

  sys.nx = 1;	   /* number of unit cells in x-direction */
  sys.ny = 1;
  sys.nz = 1;

  /* check */
  sys.N = 60;  /* total number of ions [diamond strct.] */

  sys.Lx = sys.Ax * (double)sys.nx;   /* basic MD box size */
  sys.Ly = sys.Ay * (double)sys.ny;
  sys.Lz = sys.Az * (double)sys.nz;

  ctl.kinds_of_ions = 1;   /* "C" : see also "identify_ion()" */

/*========================================================================*/

  /*--- UNIT conversion ---------------------------------------------*/
  /* sys.perMol:  [energy] -> [J][mol^{-1}]  unit conversion         */
  /* 6.0221367e5 = Avogadro constant in MD unit                      */
  /* [Mol] number of [Diamond] unit -> sys.N/3 in MD basic cell      */
  sys.perMol = 6.0221367e5/((double)(sys.N/2)); 

  /* delta_t unit conversion : 1.0 [fs] -> 1.0/4.07497263794 [second'] */
  sys.dt = delta_time_fs/4.07497263794;  
  sys.dt2 = sys.dt*sys.dt;
  sys.step_2_fsec = delta_time_fs;   /* [MD step] -> [fs] converter    */

  /* [GPa] * press_unit_conversion = [Pa']                             */
  press_unit_conversion = 1.0/(  (1.6605402*10000.0)/
                      ((4.074972637944947)*(4.074972637944947))  );

  ctl.press_X = set_press_GPa_X*press_unit_conversion;
  ctl.press_Y = set_press_GPa_Y*press_unit_conversion;
  ctl.press_Z = set_press_GPa_Z*press_unit_conversion;  

  sys.kB = 1.380658e-5;   /* Boltzmann constant [IEMD unit] */

  /* [energy] -> [temperature] = [K] */
  /* T = (2/(3NKb))Sum(E_k), where Sum(E_k)=(1/2)Sum(mv^2) */
  sys.e2t = 2.0 / (3.0 * ((double)(sys.N) * sys.kB));

  /* [Pa'] -> [GPa] */
  sys.pp2gpa = (1.6605402*10000.0)/((4.074972637944947)*(4.074972637944947));

  /*=================================================================*/
  /* TBMD specific part [Carbon]                                     */
  /*=================================================================*/
  /*-------------  H_tb ---------------------------------------------*/

  /* [eV] -> [Energy] unit converter                          */
  /* 1 [eV] = 1.60217733*10^(-19) [J] -> 0.160217733 [energy] */
  htb.eV2E = 0.160217733;

  /* energy parameters for Htb [eV] -> Htb [Energy] */
  htb.e_s = -2.99 *htb.eV2E;     /* [Energy] */
  htb.e_p = 3.71 *htb.eV2E;      /* [Energy] */
  htb.v_sssig = -5.00 *htb.eV2E; /* [Energy] */
  htb.v_spsig = 4.70 *htb.eV2E;  /* [Energy] */
  htb.v_ppsig = 5.50 *htb.eV2E;  /* [Energy] */
  htb.v_pppi = -1.55 *htb.eV2E;  /* [Energy] */

  /* parameters for function s(r) [1] */
  htb.s_n    = 2.0;      /* [dimensionless] */
  htb.s_nc   = 6.5;      /* [dimensionless] */
  htb.s_rc   = 2.18;     /* [A] */
  htb.s_r0   = 1.536329; /* [A] */
  htb.s_r1   = 2.45;     /* cutoff-1   [A]   */ 
  htb.s_r1_2 = 6.0025;   /* cutoff-1^2 [A^2] */ 
  htb.s_rm   = 2.6;      /* cutoff-2   [A]   */
  htb.s_rm_2 = 6.76;     /* cutoff-2^2 [A^2] */

  /* parameters for function s(r) [2] */
  htb.s_c0 =  6.7392620074314/1000.0; /* [dimensionless] */
  htb.s_c1 = -8.1885359517898/100.0;  /* [A^(-1)]        */
  htb.s_c2 =  0.1932365259144;        /* [A^(-2)]        */
  htb.s_c3 =  0.3542874332380;        /* [A^(-3)]        */

  /* x^y = exp(y*log(x))  [y = integer, x><0] or [y=float, x>0]  */
  htb.r0_rc_nc = exp(htb.s_nc*log(htb.s_r0/htb.s_rc));  /* =(r0/rc)^nc */

  /*-------------  Repulsion ----------------------------------------*/
  phi.phi   = 8.18555 * htb.eV2E;   /* [Energy] */
  phi.m     = 3.30304;  /* [dimensionless] */
  phi.mc    = 8.6655;   /* [dimensionless] */
  phi.dc    = 2.1052;   /* [A] */
  phi.d0    = 1.64;     /* [A] */
  phi.d1    = 2.57;     /* cutoff-1    [A]   */ 
  phi.d1_2  = 6.6049;   /* cutoff-1^2  [A^2] */ 
  phi.rm    = 2.6;      /* cutoff-2    [A]   */ 
  phi.rm_2  = 6.76;     /* cutoff-2^2  [A^2] */ 

  /* x^y = exp(y*log(x))  [y = integer, x><0] or [y=float, x>0]  */
  phi.d0_dc_mc = exp(phi.mc*log(phi.d0/phi.dc));  /* =(d0/dc)^mc */

  /* parameters for repulsive function phi(r) [2] */
  phi.c0p =  ( 2.2504290109/100000000.0);  /* [eV]       */
  phi.c1p =  (-1.4408640561/1000000.0);    /* [eV/A]     */
  phi.c2p =  ( 2.1043303374/100000.0);     /* [eV/A^2]   */
  phi.c3p =  ( 6.6024390226/100000.0);     /* [eV/A^3]   */

  phi.c0f = -2.5909765118191;              /* [eV]            */
  phi.c1f =  0.5721151498619;              /* [dimensionless] */
  phi.c2f = -1.7896349903996/1000.0;       /* [eV^(-1)]       */
  phi.c3f =  2.3539221516757/100000.0;     /* [eV^(-2)]       */
  phi.c4f = -1.24251169551587/10000000.0;  /* [eV^(-3)]       */

   /* unit conversion [eV] -> [Energy] ----------------------------------*/
   phi.c0p = phi.c0p * htb.eV2E;    /* [Energy]     */
   phi.c1p = phi.c1p * htb.eV2E;    /* [Energy/A]   */
   phi.c2p = phi.c2p * htb.eV2E;    /* [Energy/A^2] */
   phi.c3p = phi.c3p * htb.eV2E;    /* [Energy/A^3] */
  
   phi.c0f = phi.c0f * htb.eV2E;                      /* [Energy]        */
   phi.c1f = phi.c1f;                                 /* [dimensionless] */
   phi.c2f = phi.c2f / htb.eV2E;                      /* [Energy^(-1)    */
   phi.c3f = phi.c3f / (htb.eV2E*htb.eV2E);           /* [Energy^(-2)    */
   phi.c4f = phi.c4f / (htb.eV2E*htb.eV2E*htb.eV2E);  /* [Energy^(-3)    */
   /*--------------------------------------------------------------------*/

  /*=================================================================*/
}

/*****
  identify_ion:
*****/
void	identify_ion(void)
{
  short	i=0;

  /* every atom is "C" in pure Carbon system */

  for(i=0; i<sys.N; i++) { 
    sys.ion[i] = 0; 
  }
}

/*****
  void set_roc(void)
  set up the particle positions
*****/
void	set_roc(void)
{

  int i;

/****
  C_60 structure
****/

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

  for(i=0; i<60; i++){ /* for taking origin to center of cell */
    sys.rx[i] += sys.Ax/2.0;
    sys.ry[i] += sys.Ay/2.0;
    sys.rz[i] += sys.Az/2.0;
  }
}

void set_potential(void)
{
  double m_C;  /* mass of Carbon atom */

  m_C = 12.011 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */
  ion.m[0] = m_C;

}

