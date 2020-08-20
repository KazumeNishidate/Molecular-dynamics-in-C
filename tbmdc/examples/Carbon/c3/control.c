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

  ctl.calc_max = 10000;  /* maximum MD time step                     */
  delta_time_fs = 0.3;   /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 300.0;       /* Temperature setting [K] */
  ctl.t_control_step = 5;
  set_press_GPa_X = 0.01;  /* Pressure setting  [GPa] */
  set_press_GPa_Y = 0.01; 
  set_press_GPa_Z = 0.01; 
  ctl.p_control_step = 20;
  sys.pres_scalling_factor = 0.5;

  /* check */
  sys.Ax = 1.31;     /* initial lattice constant */
  sys.Ay = 1.31;     /* initial lattice constant */
  sys.Az = 1.31;     /* initial lattice constant */
/*
  sys.Ay = 3.55984;
  sys.Az = 3.55984;
*/
  sys.nx = 1;	   /* number of unit cells in x-direction */
  sys.ny = 1;
  sys.nz = 1;

  /* check */
  sys.N = 3*sys.nx*sys.ny*sys.nz;  /* total number of ions [diamond strct.] */

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

  /* tight-binding matrix memory allocation */

  /* Htb */
  htb.mat1  = (double *)calloc(sys.N*sys.N*4*4, sizeof(double));

  /* FH force matrix */
  htb.mat2x = (double *)calloc(sys.N*sys.N*4*4, sizeof(double));
  htb.mat2y = (double *)calloc(sys.N*sys.N*4*4, sizeof(double));
  htb.mat2z = (double *)calloc(sys.N*sys.N*4*4, sizeof(double));

  /* store the occupied electron IDs */
  htb.occupied = (int *)calloc(sys.N*4, sizeof(int));

  /* working valiables to sort the eigen values in order */
  wv.eigen_values = (double *)calloc(sys.N*4, sizeof(double));
  wv.eigen_values_address = (int *)calloc(sys.N*4, sizeof(int));

  /* Jacobi unitary transformation --------------------------------*/

   /* unitary matrix */
   jacobi.Unitary = (double *)calloc(sys.N*sys.N*16, sizeof(double));

   /* U^(T) * Htb */
   jacobi.Ut_tb = (double *)calloc(sys.N*sys.N*16, sizeof(double));

   /* U^(T) * F * U */
   jacobi.fx = (double *)calloc(sys.N*sys.N*16, sizeof(double));
   jacobi.fy = (double *)calloc(sys.N*sys.N*16, sizeof(double));
   jacobi.fz = (double *)calloc(sys.N*sys.N*16, sizeof(double));

  /*---------------------------------------------------------------*/

  /* lookup list of repulsive potential for each atom */
  phi.sum_phi        = (double *)calloc(sys.N, sizeof(double));
  phi.f_dash_sum_phi = (double *)calloc(sys.N, sizeof(double));

  phi.phi_dash_x     = (double *)calloc(sys.N*sys.N, sizeof(double));
  phi.phi_dash_y     = (double *)calloc(sys.N*sys.N, sizeof(double));
  phi.phi_dash_z     = (double *)calloc(sys.N*sys.N, sizeof(double));

  phi.sum_phi_dash_x     = (double *)calloc(sys.N, sizeof(double));
  phi.sum_phi_dash_y     = (double *)calloc(sys.N, sizeof(double));
  phi.sum_phi_dash_z     = (double *)calloc(sys.N, sizeof(double));

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
  int i=0;

  /* 3-Carbon system: bond length at equilibrium is 1.310 [A]           */
  /* total number of particles sys.N was set in "get_control_param()"   */

  /* MD-box size */
  sys.Lx = 14.0;
  sys.Ly = 12.0;
  sys.Lz = 12.0;

  sys.rx[i] = 6.0;
  sys.ry[i] = 6.0;
  sys.rz[i] = 6.0;
  i++;

  sys.rx[i] = 6.0+sys.Ax;
  sys.ry[i] = 6.0;
  sys.rz[i] = 6.0;
  i++;

  sys.rx[i] = 6.0+sys.Ax+sys.Ax;
  sys.ry[i] = 6.0;
  sys.rz[i] = 6.0;
  i++;

}

void set_potential(void)
{
  double m_C;  /* mass of Carbon atom */

  m_C = 12.011 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */
  ion.m[0] = m_C;

}

