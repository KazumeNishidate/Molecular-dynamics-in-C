#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

void   init_mem_chebyshev_coeff(void){

  ch.Npl = 101;  /* take the polynomials up to 50th  order */
  ch.Emax = 40.0*htb.eV2E;

  fermi.electrons = (double)sys.N*2.0; /* 2 electrons per 1 Carbon atom */

  /* initial chemical potentail setting */
  ch.myu = 3.71* htb.eV2E;

  /*-----------------------------------------------------------------*/
  /* localized region for the target atom */
  local.Rc = 5.7; /* localization cut off radius for Htb */
  local.Rc2 = local.Rc*local.Rc;

  local.maxN = 63;   /* possible max number of atoms in local region */

  /* store the NN atoms in R_cutoff */
  local.lookup = (int *)calloc(sys.N*local.maxN+1, sizeof(int));
  local.address = (int *)calloc(sys.N+1, sizeof(int));

  /* total number of neighbor atoms for each atom */
  local.num_of_nn_atoms = (int *)calloc(sys.N+1, sizeof(int));

  /* ID list of atoms to make up the LOCAL Htb at hamiltonian.c */
  local.atoms = (int *)calloc(local.maxN+1, sizeof(int));
  /*-----------------------------------------------------------------*/

  ch.H  = (double *)calloc(ch.Npl+2, sizeof(double));
  ch.F  = (double *)calloc(ch.Npl+2, sizeof(double));
  fermi.T1 = (double *)calloc(local.maxN*local.maxN*4*4+1, sizeof(double));
  fermi.T2 = (double *)calloc(local.maxN*local.maxN*4*4+1, sizeof(double));
  fermi.T3 = (double *)calloc(local.maxN*local.maxN*4*4+1, sizeof(double));
  fermi.htb = (double *)calloc(local.maxN*local.maxN*4*4+1, sizeof(double));
  fermi.mat = (double *)calloc(local.maxN*local.maxN*4*4+1, sizeof(double));

  /*--------------------------------------------------------------------*/
  /* tight-binding matrix memory allocation */

  /* Htb */
  htb.mat1  = (double *)calloc(local.maxN*local.maxN*4*4, sizeof(double));

  /* FH force matrix */
  htb.mat2x = (double *)calloc(local.maxN*local.maxN*4*4, sizeof(double));
  htb.mat2y = (double *)calloc(local.maxN*local.maxN*4*4, sizeof(double));
  htb.mat2z = (double *)calloc(local.maxN*local.maxN*4*4, sizeof(double));

  /* store the occupied electron IDs */  
  /* checklocal [needless?] */
  htb.occupied = (int *)calloc(local.maxN*2, sizeof(int));

  htb.potential = (double *)calloc(sys.N, sizeof(double));

  /*---------------------------------------------------------------*/

  /* lookup list of repulsive potential for each atom */

  /*---------------------------------------------------------------*/

  /* lookup list of repulsive potential for each atom */
  phi.sum_phi        = (double *)calloc(sys.N, sizeof(double));
  phi.f_dash_sum_phi = (double *)calloc(sys.N, sizeof(double));

  phi.phi_dash_x     = (double *)calloc(sys.N*sys.N, sizeof(double));
  phi.phi_dash_y     = (double *)calloc(sys.N*sys.N, sizeof(double));
  phi.phi_dash_z     = (double *)calloc(sys.N*sys.N, sizeof(double));

  /*=================================================================*/

  /* ----- allocate memory to member of PRESSURE structure ----- */
  /* -----         taka (Wednesday, June 7, 2000)          ----- */
  press.dx = (double *)calloc(sys.N*sys.N+1,sizeof(double));
  press.dy = (double *)calloc(sys.N*sys.N+1,sizeof(double));
  press.dz = (double *)calloc(sys.N*sys.N+1,sizeof(double));

  press.matX = (double *)calloc(local.maxN*local.maxN*4*4, sizeof(double));
  press.matY = (double *)calloc(local.maxN*local.maxN*4*4, sizeof(double));
  press.matZ = (double *)calloc(local.maxN*local.maxN*4*4, sizeof(double));

}

