#include	<stdio.h>
#include	<stdlib.h> 
#include	<math.h>

#define PI 3.141592654
#define PI2 6.283185308

FILE *fpna, *fpcl, *fp_positions, *fpout,*frho_na,*frho_s_c;

int num_of_particles;     /*  total number of particles    */
int num_of_ion0;
int num_of_ion1;
int total_time_step;      /*  number of maximum time step  */
double delta_t;           /*  time increment delta_t [fs]  */
double hh, kk, ll;        /*  reciprocal lattice point     */

double Z_ion0, Z_ion1;    /*  weighting factor [charge]                 */
double B_ion0, B_ion1;    /*  weighting factor [scattering length]      */
double M_ion0, M_ion1;    /*  weighting factor [mass]                   */

double delta_w;
double *gskw;
double *skw;
double *skw_PP, *skw_PN, *skw_NN;
double *rho_k_c_Na, *rho_k_s_Na;
double *rho_k_c_Cl, *rho_k_s_Cl;
double *rho_k_c_Na_rft, *rho_k_s_Na_rft;
double *rho_k_c_Cl_PN, *rho_k_s_Cl_PN;
double *rho_k_c_Na_Im, *rho_k_s_Na_Im;
double *rho_k_c_Cl_Im, *rho_k_s_Cl_Im;


