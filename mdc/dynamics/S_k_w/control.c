#include  "sk.h"

void  init_param(void)
{
  double Ax, Ay, Az;   /* unit cell size */
  double rh, rk, rl;   /* reciprocal lattice point {h, k, l}  */

/*-----------------------------------------------------------------*/

  Ax = 5.63;   /* lattice constant in [A] unit: [NaCl] */
  Ay = 5.63;
  Az = 5.63;

  num_of_particles = 216;   /*  total number of particles      */
  num_of_ion0 = 108;        /*  PLUS charge */     
  num_of_ion1 = 108;        /*  MINUS charge */

  total_time_step = 8192;    /*  number of maximum MD time step */
  delta_t = 1.0;            /*  delta-t [fs] */

  rh = 1.0/3.0;  /* reciprocal lattice point {h} in [2PI/Ax unit] */
  rk = 6.0/3.0;      /* reciprocal lattice point {k} in [2PI/Ay unit] */
  rl = 0.0;      /* reciprocal lattice point {l} in [2PI/Az unit] */

  Z_ion0 =  1.0;   /* effective charge of Na ion [arbitrary unit]  */
  Z_ion1 = -1.0;   /* effective charge of Cl ion [arbitrary unit]  */
  B_ion0 =  0.52;  /* scattering length of Na ion [arbitrary unit] */
  B_ion1 =  1.47;  /* scattering length of Cl ion [arbitrary unit] */
  M_ion0 =  23.0;  /* mass of Na ion */
  M_ion1 =  35.5;  /* mass of Cl ion */

/*-----------------------------------------------------------------*/

  hh = rh * 2.0 * PI/Ax;
  kk = rk * 2.0 * PI/Ay;
  ll = rl * 2.0 * PI/Az;

  /* delta-omega  [rad]/[ps] */
  delta_w = (1000.0*PI2)/((double)total_time_step*delta_t);  

}

void  init_mem(void) /*--- dynamic memory allocation ---*/
{
  rho_k_c_Na = (double *)calloc(total_time_step, sizeof(double));
  rho_k_c_Cl = (double *)calloc(total_time_step, sizeof(double));
  rho_k_s_Na = (double *)calloc(total_time_step, sizeof(double));
  rho_k_s_Cl = (double *)calloc(total_time_step, sizeof(double));
  rho_k_c_Na_rft = (double *)calloc(total_time_step, sizeof(double));
  rho_k_s_Na_rft = (double *)calloc(total_time_step, sizeof(double));
  rho_k_c_Cl_PN = (double *)calloc(total_time_step, sizeof(double));
  rho_k_s_Cl_PN = (double *)calloc(total_time_step, sizeof(double));
  rho_k_c_Na_Im = (double *)calloc(total_time_step, sizeof(double));
  rho_k_s_Na_Im= (double *)calloc(total_time_step, sizeof(double));
  rho_k_c_Cl_Im = (double *)calloc(total_time_step, sizeof(double));
  rho_k_s_Cl_Im = (double *)calloc(total_time_step, sizeof(double));

  skw_PP = (double *)calloc(total_time_step, sizeof(double));
  skw_NN = (double *)calloc(total_time_step, sizeof(double));
  skw_PN = (double *)calloc(total_time_step, sizeof(double));
  skw = (double *)calloc(total_time_step, sizeof(double));
  gskw = (double *)calloc(total_time_step, sizeof(double));
}

