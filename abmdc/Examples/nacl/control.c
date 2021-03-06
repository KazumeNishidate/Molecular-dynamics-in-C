#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

/*************************************************************/
/*  [ setting note ]                                         */
/*  NaCl crystal system                                      */
/*  one unit cell with 4 Na and 4Cl atoms                    */
/*                                 2000. Jan. project IEMD   */
/*************************************************************/

void  set_param(void)
{
  wv.E_cut  = 7.0;   /* plane wave cut off energy [Hartree]               */
  wv.E_cuti = 2.0;    /* cutoff energy at initial stage of diagonalization */
                      /* E_cut > E_cuti                                    */

  /*  wv.E_cut  = 15.0;  wv.E_cuti = 3.0;  */

  /*--- approach to the Born-Oppenheimer potential surface ---------------*/
  cpq.dt = 4.0;            /* = delta t [a.u.]                */
  cpq.myu = 200.0;	   /* pseudo mass of electron [a.u.]  */
  cpq.residual = 1.0e-5;   /* allowed residual                */
  cpq.max_qstep = 1000;    /* MAX quench step                 */
  cpq.quench = 0.5;	   /* quench rate for wave functios   */
  cpq.qstep_interval = 1;  /* interval for the quenching      */
  cpq.dt2myu = cpq.dt*cpq.dt/cpq.myu;

  /*--- used at CP MD loop -----------------------------------------------*/
  cpa.dt = 4.0;            /* = delta t [a.u.]                */
  cpa.myu = 200.0;	   /* pseudo mass [a.u.]              */
  cpa.residual = 1.0e-4;   /* allowed residual                */
  cpa.max_step = 5000;	   /* MAX MD step                     */
  cpa.quench = 1.0;	   /* quench rate for atoms motion    */
  cpa.qstep_interval = 1;  /* interval for the quenching      */
  cpa.dt2myu = cpa.dt*cpa.dt/cpa.myu;
  /*----------------------------------------------------------------------*/

  /*--- control temperature ----------------------------------------------*/
  ctl.temp = 1200.0;
  ctl.tstep = 10;
  /*----------------------------------------------------------------------*/

  ft.nx = FFTPTSX; /* FFT mesh number : defined in "heders/size.h" */
  ft.ny = FFTPTSY;
  ft.nz = FFTPTSZ;
  ft.mesh = ft.nx*ft.ny*ft.nz;
}

/*****
 *  set_atoms:
 *****/
void	set_atoms(void)
{
  int i, type;
  double z=0.0;

  atom.ntypes = 2;  /* total types of atoms = 2 */

  type = 0;  /*--- [Na] ----------------------------------------------------*/
  atom.m[type]        = 22.989768/(AVOGAD*ME);  /* atomic mass of Na [au]   */
  atom.atomic_n[type] = 11;                   /* atomic number Na = 11      */
  atom.Zv[type]       = 1.0;                  /* valence charge Zv = 1      */
  atom.nvband[type]   = (int)atom.Zv[type];   /* number of valence band     */
                             /*-- check : atom.nvband != atom.Zv in general */
  atom.norbit[type]   = 3;                    /* number of AO : L = 0, 1, 2 */
  atom.ns_orbit[type] = 0;                    /* spin-orbital = none [Na]   */
  /*------------------------------------------------------------------------*/

  type = 1;  /*--- [Cl] ----------------------------------------------------*/
  atom.m[type]        = 35.4527/(AVOGAD*ME);  /* atomic mass of Cl [au]     */
  atom.atomic_n[type] = 17;                    /* atomic number Cl = 1      */
  atom.Zv[type]       = 7.0;                  /* valence charge Zv = 7      */
  atom.nvband[type]   = (int)atom.Zv[type];   /* number of valence band     */
                             /*-- check : atom.nvband != atom.Zv in general */
  atom.norbit[type]   = 3;                    /* number of AO : L = 0, 1, 2 */
  atom.ns_orbit[type] = 0;  /* not used */    /* spin-orbital = none [Cl]   */
  /*------------------------------------------------------------------------*/

  set_rocs();       /* set up the atomic positions */
  init_velocity();  /* set up the atomic velocity  */

  for(i=0;i<atom.N;i++) { /* evaluate the total number of bands */
    z += (double)atom.nvband[atom.type[i]];
  }
  atom.nvband_all = ceil(z/2.0); /* neglecting the SPIN AO */
  /* total no. of bands = (int) (total valence charge)/2.0 */

}

/*****
 *  set_rocs:
 *  set up the positions of atoms
 *****/
void	set_rocs(void)
{
  int i;

  cell.Ax = 5.63*au.r;  /* lattice constant [NaCl rock salt] */
  cell.Ay = 5.63*au.r;  /* [angstrom] -> [atomic unit]       */          
  cell.Az = 5.63*au.r;

  cell.nx = 1;          /* number of unit cells in each direction */
  cell.ny = 1;
  cell.nz = 1;

  cell.Lx = cell.Ax * (double)cell.nx;   /* basic MD box size */
  cell.Ly = cell.Ay * (double)cell.ny;
  cell.Lz = cell.Az * (double)cell.nz;

  cell.vol = cell.Lx*cell.Ly*cell.Lz;    /* volume of the MD box */ 

  /*=====================================================================*/
  atom.N = 8;    /* total number of atoms [NaCl] */

  /* [Na] */
  cell.rx[0] = 0.5;   cell.rx[1] = 0.0;  cell.rx[2] = 0.0;   cell.rx[3] = 0.5;
  cell.ry[0] = 0.5;   cell.ry[1] = 0.0;  cell.ry[2] = 0.5;   cell.ry[3] = 0.0;
  cell.rz[0] = 0.5;   cell.rz[1] = 0.5;  cell.rz[2] = 0.0;   cell.rz[3] = 0.0;
  atom.type[0] = 0;   atom.type[1] = 0;  atom.type[2] = 0;   atom.type[3] = 0;

  /* [Cl] */
  cell.rx[4] = 0.0;   cell.rx[5] = 0.5;  cell.rx[6] = 0.5;   cell.rx[7] = 0.0;
  cell.ry[4] = 0.0;   cell.ry[5] = 0.5;  cell.ry[6] = 0.0;   cell.ry[7] = 0.5;
  cell.rz[4] = 0.0;   cell.rz[5] = 0.0;  cell.rz[6] = 0.5;   cell.rz[7] = 0.5;
  atom.type[4] = 1;   atom.type[5] = 1;  atom.type[6] = 1;   atom.type[7] = 1;

  for(i=0;i<8;i++) {
    cell.rx[i] *= cell.Lx;
    cell.ry[i] *= cell.Ly;
    cell.rz[i] *= cell.Lz;
  }

  /*=====================================================================*/
  
}

/*****
 *  initialize the BHS pseudo potential:
 *  - Ref: Phys. Rev. B., Vol 26, No. 8, pp. 4199-4228,
 *    G.B. Bachelet, D.R. Hamann, and M. Schlu"ter,
 *    "Pseudopotentials that work: From H to Pu"
 *  - Errata: Phys. Rev. B., Vol 29, No. 4, pp. 2309
 *    G.B. Bachelet, D.R. Hamann, and M. Schlu"ter,
*****/
void   init_bhs_pseudo_potential(void)
{
  int type, eru;

  /*--- Na BHS pseudo potentail parameters ------------------------------*/

  type = 0; /* Na */

  /* core part */
  bhs.core_orbit_alpha[type][0] = 1.71;  /* alpha_core_1 */
  bhs.core_orbit_alpha[type][1] = 0.50;  /* alpha_core_2 */
  bhs.core_orbit_c[type][0] = 5.1815;    /* c_core_1     */
  bhs.core_orbit_c[type][1] = -4.1815;   /* c_core_2     */

  /* orbital part */
  eru = 0;
  bhs.alpha[type][eru][1] = 0.99;
  bhs.alpha[type][eru][2] = 1.10;
  bhs.alpha[type][eru][3] = 1.24;
  bhs.c[type][eru][0] = -2.4718;  bhs.c[type][eru][1] = 0.3334;
  bhs.c[type][eru][2] = 0.0619;   bhs.c[type][eru][3] = 0.0890;
  bhs.c[type][eru][4] = -0.0014;  bhs.c[type][eru][5] = -0.0123;

  eru = 1;
  bhs.alpha[type][eru][1] = 0.51;
  bhs.alpha[type][eru][2] = 0.65;
  bhs.alpha[type][eru][3] = 0.84;
  bhs.c[type][eru][0] = -1.6202;  bhs.c[type][eru][1] = -0.4908;
  bhs.c[type][eru][2] = -0.0861;  bhs.c[type][eru][3] = 0.0375;
  bhs.c[type][eru][4] = -0.0161;  bhs.c[type][eru][5] = 0.0070;

  eru = 2;
  bhs.alpha[type][eru][1] = 0.38;
  bhs.alpha[type][eru][2] = 0.55;
  bhs.alpha[type][eru][3] = 0.73;
  bhs.c[type][eru][0] = -0.9415;  bhs.c[type][eru][1] = -0.9710;
  bhs.c[type][eru][2] = -0.2336;  bhs.c[type][eru][3] = -0.0593;
  bhs.c[type][eru][4] = -0.0228;  bhs.c[type][eru][5] = 0.0455;
  /*---------------------------------------------------------------------*/

  /*--- Cl BHS pseudo potentail parameters ------------------------------*/

  type = 1; /* Cl */

  /* core part */
  bhs.core_orbit_alpha[type][0] = 3.48;  /* alpha_core_1 */
  bhs.core_orbit_alpha[type][1] = 1.38;  /* alpha_core_2 */
  bhs.core_orbit_c[type][0] = 1.3860;    /* c_core_1     */
  bhs.core_orbit_c[type][1] = -0.3860;   /* c_core_2     */

  /* orbital part */
  eru = 0;
  bhs.alpha[type][eru][1] = 4.94;
  bhs.alpha[type][eru][2] = 9.61;
  bhs.alpha[type][eru][3] = 15.05;
  bhs.c[type][eru][0] = -3.6651;  bhs.c[type][eru][1] = 1.2609;
  bhs.c[type][eru][2] = -0.5528;  bhs.c[type][eru][3] = -0.3237;
  bhs.c[type][eru][4] = 0.0368;  bhs.c[type][eru][5] = 0.0150;

  eru = 1;
  bhs.alpha[type][eru][1] = 2.41;
  bhs.alpha[type][eru][2] = 3.16;
  bhs.alpha[type][eru][3] = 4.73;
  bhs.c[type][eru][0] = -2.3089;  bhs.c[type][eru][1] = -0.0556;
  bhs.c[type][eru][2] = 0.0784;  bhs.c[type][eru][3] = 0.0357;
  bhs.c[type][eru][4] = -0.0080;  bhs.c[type][eru][5] = -0.0060;

  eru = 2;
  bhs.alpha[type][eru][1] = 4.04;
  bhs.alpha[type][eru][2] = 4.83;
  bhs.alpha[type][eru][3] = 5.40;
  bhs.c[type][eru][0] = 0.0968;  bhs.c[type][eru][1] = -0.6838;
  bhs.c[type][eru][2] = -0.1482;  bhs.c[type][eru][3] = -0.3090;
  bhs.c[type][eru][4] = 0.0218;  bhs.c[type][eru][5] = 0.0735;
  /*---------------------------------------------------------------------*/

  for(type=0;type<atom.ntypes;type++) {
    init_bhs_c2a(type); /* evaluate the Ai cofficients */
  }

}

