#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_ 1

/*********************************************************************/
/* RUN OPTION setting                                                */
/*********************************************************************/
/* Chebyshev polynomials expansion of Fermi matrix                   */
#define chebyshev
/*********************************************************************/


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

  double *potential;

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
  double *phi_dash_x;
  double *phi_dash_y;
  double *phi_dash_z;

} repulsion;

typedef struct{  /* coefficients for the Chebyshev polynomials */
  int    Npl; /* number of expansion terms */
  double myu; /* chemical potential */ 
  double *H;  /* chebyshev coefficients for Htb expansion */
  double *F;  /* chebyshev coefficients for Ftb expansion */
  double Emax;

} chebyshev_coeff;

typedef struct{  /* Fermi matrix */

  double *T1; /* Fermi matrix to use in Chebyshev polynomials */
  double *T2;
  double *T3;
  double *htb; /* htb = Htb/Emax */
  double *mat;
  double number_of_electrons;  /* = Tr[Fermi(htb)] */
  double electrons;

} fermi_matrix;

typedef struct{

  int maxN;      /* possible max number of atoms in LOCAL region */

  double Rc;     /* cut-off distance of the Local Hamiltonian */
  double Rc2;

  int *lookup;   /* nn atom lookup table list */
  int *address;  /* store the address of the local.lookup */

  int *num_of_nn_atoms; /* total number of neighbor atoms for each atom */

  int *atoms;  /* ID list of atoms to make up the LOCAL Htb */
  
} Order_N;

typedef struct tagPRESSURE{
  double *dx, *dy, *dz;
  double *matX, *matY, *matZ;
  double virX, virY, virZ;
}PRESSURE;

typedef struct{

  int    ch_n;        /* <--- index number of chiral vector */
  int    ch_m;        /* <--- the same as above */
  int    *chain;      /* <--- number of atoms in one "chain" */
  int    atoms;       /* <--- number of atoms in one "tube" */
  double norm_ch;     /* <--- norm of chiral vector */
  double bond_length; /* <--- bond-length between carbon atoms */
  double ch_theta;    /* <--- chiral angle */
  double rad;         /* <--- radius of tube */
  double T;           /* <--- norm of basic translational vector */
  double *rx, *ry, *rz;

} NANOTUBE;

typedef struct tagMULTIWALL{
    int num_of_wall;
    int atoms;
    double T;
    double rad;
    double *rx, *ry, *rz;
} MULTIWALL;


/*------------------- declaration for the structures ----------------------*/
tight_binding_hamiltonian htb;
repulsion phi;
chebyshev_coeff ch;
fermi_matrix fermi;
Order_N local;

PRESSURE press;

NANOTUBE tube;
MULTIWALL multi;

#endif
