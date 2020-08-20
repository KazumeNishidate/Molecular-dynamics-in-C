/*********************************************************************/
/*  Tight Binding Hamiltonian part                                   */
/*********************************************************************/

typedef struct{

  double e_s;    /* energy parameters for Htb */
  double e_p;
  double v_sssig;
  double v_spsig;
  double v_ppsig;
  double v_pppi;

  double s_n;    /* parameters for function s(r) [1] */
  double s_nc;
  double s_rc;   /* [A] */
  double s_r0;   /* [A] */
  double s_r1;   /* [A] cutoff-1   */ 
  double s_r1_2; /* [A] cutoff-1^2 */ 
  double s_rm;   /* [A] cutoff-2   */
  double s_rm_2; /* [A] cutoff-2^2 */

  double s_c0;   /* parameters for function s(r) [2] */
  double s_c1;
  double s_c2;
  double s_c3;
  double r0_rc_nc;  /* =(r0/rc)^nc */

  double *mat1;   /* matrix for the potential */
  double *mat2x;  /* matrix for the force in x-direction */
  double *mat2y;  /* matrix for the force in y-direction */
  double *mat2z;  /* matrix for the force in z-direction */

  double eV2E;   /* 1 [eV] -> 0.160217733 [energy] */

  int *occupied; /* address for occupied 2 electrons for each atom */

} tight_binding_hamiltonian;

typedef struct{
  double phi;   /* parameters for repulsive function phi(r) [1] */
  double m;
  double mc;
  double dc;
  double d0;
  double d1;
  double d1_2;
  double rm;
  double rm_2;
  double d0_dc_mc;  /* =(d0/dc)^mc */

  double c0p;   /* parameters for repulsive function phi(r) [2] */
  double c1p;
  double c2p;
  double c3p;

  double c0f;
  double c1f;
  double c2f;
  double c3f;
  double c4f;

  double *sum_phi;
  double *f_dash_sum_phi;
  double *sum_phi_dash_x;
  double *sum_phi_dash_y;
  double *sum_phi_dash_z;
  double *phi_dash_x;
  double *phi_dash_y;
  double *phi_dash_z;

} repulsion;

typedef struct{

  double *eigen_values; /* eigen values list to sort */
  int    *eigen_values_address;

} working_valiables;

typedef struct{
  double *Unitary;

  double *Ut_tb;

  double *fx;
  double *fy;
  double *fz;
} jacobi_working_valiables;


/* ----- taka (Monday, June 5, 2000) ----- */
typedef struct tagPRESSURE{
  double *dx, *dy, *dz;
  double *matX, *matY, *matZ;
  double virX, virY, virZ;
}PRESSURE;



/*------------------- declaration for the structures ----------------------*/
tight_binding_hamiltonian htb;
repulsion phi;
working_valiables wv;
jacobi_working_valiables jacobi;

PRESSURE press; /* ----- taka (Monday, June 5, 2000) ----- */

