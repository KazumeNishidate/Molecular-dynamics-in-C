#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

/*************************************************************/
/*  [ setting note ]                                         */
/*  V2O5                                                     */
/*                                 1999. Dec. project IEMD   */
/*************************************************************/

void  set_param(void)
{
  wv.E_cut  = 4.0;   /* plane wave cut off energy [Hartree]               */
  wv.E_cuti = 1.5;   /* cutoff energy at initial stage of diagonalization */
                     /* E_cut > E_cuti                                    */

  /*--- approach to the Born-Oppenheimer potential surface ---------------*/
  cpq.dt = 2.0;            /* = delta t [a.u.]                */
  cpq.myu = 50.0;	   /* pseudo mass of electron [a.u.]  */
  cpq.residual = 1.0e-3;   /* allowed residual                */
  cpq.max_qstep = 10000;   /* MAX quench step                 */
  cpq.quench = 0.5;	   /* quench rate for wave functios   */
  cpq.qstep_interval = 1;  /* interval for the quenching      */
  cpq.dt2myu = cpq.dt*cpq.dt/cpq.myu;

  /*--- used at CP MD loop -----------------------------------------------*/
  cpa.dt = 2.0;            /* = delta t [a.u.]                */
  cpa.myu = 100.0;	   /* pseudo mass [a.u.]              */
  cpa.residual = 5.0e-3;   /* allowed residual                */
  cpa.max_step = 30000;	   /* MAX MD step                     */
  cpa.quench = 1.0;	   /* quench rate for atoms motion    */
  cpa.qstep_interval = 1;  /* interval for the quenching      */
  cpa.dt2myu = cpa.dt*cpa.dt/cpa.myu;
  /*----------------------------------------------------------------------*/

  /*--- control temperature ----------------------------------------------*/
  ctl.temp = 10.0;
  ctl.step = 5;
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
  int i;
  double z=0.0;

  set_rocs();       /* set up the atomic positions */
  init_velocity();  /* set up the atomic velocity  */

  atom.ntypes = 2;  /* total types of atoms = 1 (monoatomic system) */

  /* for type V atom [V] ------------------------------------------------*/
  atom.m[0]        = 50.9415/(AVOGAD*ME);  /* atomic mass of V [au]      */
  atom.atomic_n[0] = 23;                   /* atomic number V = 23       */
  atom.Zv[0]       = 5.0;                  /* valence charge Zv = 5      */
  atom.nvband[0]   = (int)atom.Zv[0];      /* number of valence band     */
                          /*-- check : atom.nvband != atom.Zv in general */
  atom.norbit[0]   = 3;                    /* number of AO : L = 0, 1, 2 */
  atom.ns_orbit[0] = 2;                    /* spin-orbital = 2           */
  /*---------------------------------------------------------------------*/

  /* for type 1 atom [O] ------------------------------------------------*/
  atom.m[1]        = 15.9994/(AVOGAD*ME);  /* atomic mass of H [au]      */
  atom.atomic_n[1] = 8;                    /* atomic number H = 1        */
  atom.Zv[1]       = 6.0;                  /* valence charge Zv = 1      */
  atom.nvband[1]   = (int)fabs((atom.Zv[1])); /* number of valence band     */
                          /*-- check : atom.nvband != atom.Zv in general */
  atom.norbit[1]   = 3;                    /* number of AO : L = 0, 1, 2 */
  atom.ns_orbit[1] = 0;                    /* spin-orbital = none [O]    */
  /*---------------------------------------------------------------------*/

  for(i=0;i<4;i++){
    atom.type[i] = 0;
  }
  for(i=4;i<14;i++){
    atom.type[i] = 1;
  }

  for(i=0;i<atom.N;i++) { /* evaluate the total number of bands */
    z += (double)atom.nvband[atom.type[i]];
  }
  atom.nvband_all = ceil(z/2.0); /* neglecting the SPIN AO  */
  /* total no. of bands = (int) (total valence charge)/2.0 */

}

/*****
 *  set_rocs:
 *  set up the positions of atoms
 *****/
void	set_rocs(void)
{
  int i;

  cell.Ax = 11.48*au.r;  /* lattice constant [Si diamond]     */
  cell.Ay = 4.36*au.r;   /* [angstrom] -> [atomic unit]       */
  cell.Az = 3.55*au.r;

  cell.nx = 1;      /* number of unit cells in each direction */
  cell.ny = 1;
  cell.nz = 1;

  cell.Lx = cell.Ax * (double)cell.nx;   /* basic MD box size */
  cell.Ly = cell.Ay * (double)cell.ny;
  cell.Lz = cell.Az * (double)cell.nz;

  cell.vol = cell.Lx*cell.Ly*cell.Lz;    /* volume of the MD box */ 

  atom.N = 14;

  /* normalized coordinates of atoms */
  /* V */
  cell.rx[0] = 0.855;   cell.rx[2] = 6.615;
  cell.ry[0] = 2.639;   cell.ry[2] = 1.721;
  cell.rz[0] = 0.887;   cell.rz[2] = 2.667;
		                           
  cell.rx[1] = 3.185;   cell.rx[3] = 8.945;
  cell.ry[1] = 1.721;   cell.ry[3] = 2.639;                    
  cell.rz[1] = 2.667;   cell.rz[3] = 0.887;

  /* O */
  cell.rx[4] = 0.855;   cell.rx[9] = 6.975; 
  cell.ry[4] = 4.180;   cell.ry[9] = 2.180; 
  cell.rz[4] = 0.887;   cell.rz[9] = 0.887; 
		                            
  cell.rx[5] = 3.185;   cell.rx[10] = 1.215;
  cell.ry[5] = 0.180;   cell.ry[10] = 2.180;
  cell.rz[5] = 2.667;   cell.rz[10] = 2.667;
		                            
  cell.rx[6] = 6.615;   cell.rx[11] = 8.585;
  cell.ry[6] = 0.180;   cell.ry[11] = 2.180;
  cell.rz[6] = 2.667;   cell.rz[11] = 2.667;
		                            
  cell.rx[7] = 8.945;   cell.rx[12] = 10.655;
  cell.ry[7] = 4.180;   cell.ry[12] = 2.180;
  cell.rz[7] = 0.887;   cell.rz[12] = 0.887;
		                            
  cell.rx[8] = 2.825;   cell.rx[13] = 4.895;
  cell.ry[8] = 2.180;   cell.ry[13] = 2.180;
  cell.rz[8] = 0.887;   cell.rz[13] = 2.667;

  for(i=0; i<atom.N; i++) {
    cell.rx[i] *= au.r;
    cell.ry[i] *= au.r;
    cell.rz[i] *= au.r;
  }
  
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

  /*--- V BHS pseudo potentail parameters ------------------------------*/

  type = 0; /* V monoatomic system */

  /* core part */
  bhs.core_orbit_alpha[type][0] = 5.14;  /* alpha_core_1 */
  bhs.core_orbit_alpha[type][1] = 1.11;  /* alpha_core_2 */
  bhs.core_orbit_c[type][0] = 2.9680;    /* c_core_1     */
  bhs.core_orbit_c[type][1] = -1.9680;   /* c_core_2     */

  /* orbital part */
  eru = 0;
  bhs.alpha[type][eru][1] = 1.23;
  bhs.alpha[type][eru][2] = 1.54;
  bhs.alpha[type][eru][3] = 1.97;
  bhs.c[type][eru][0] = -6.6485;
  bhs.c[type][eru][1] = -0.3951; 
  bhs.c[type][eru][2] = -0.5795;
  bhs.c[type][eru][3] =  0.1797; 
  bhs.c[type][eru][4] = -0.0151;
  bhs.c[type][eru][5] =  0.0176;

  eru = 1;
  bhs.alpha[type][eru][1] = 1.03;
  bhs.alpha[type][eru][2] = 1.33;
  bhs.alpha[type][eru][3] = 1.43;
  bhs.c[type][eru][0] = -5.5539;
  bhs.c[type][eru][1] = -0.9047;
  bhs.c[type][eru][2] = -0.8758;
  bhs.c[type][eru][3] =  0.1629;
  bhs.c[type][eru][4] = -0.1369;
  bhs.c[type][eru][5] =  0.0810;

  eru = 2;
  bhs.alpha[type][eru][1] = 4.60;
  bhs.alpha[type][eru][2] = 3.28;
  bhs.alpha[type][eru][3] = 17.66;
  bhs.c[type][eru][0] = 2.4040;
  bhs.c[type][eru][1] = -1.2280;
  bhs.c[type][eru][2] = 0.4141;
  bhs.c[type][eru][3] = -0.4851;
  bhs.c[type][eru][4] = -0.2405;
  bhs.c[type][eru][5] = 0.3451;

  /* SO-1 */
  /* 0.57
     1.30
     1.36
     -0.0221
     0.0064
     -0.0025
     -0.0011
     -0.0001
     0.0000
     */
  /* SO-2 */
  /* 9.12
     12.58
     16.56
     -0.0074
     -0.0022
     0.0011
     -0.0004
     0.0000
     0.0002
     */

  /*---------------------------------------------------------------------*/

  /*--- O BHS pseudo potentail parameters -------------------------------*/

  type = 1; /* O */

  /* core part */
  bhs.core_orbit_alpha[type][0] = 18.09;  /* alpha_core_1 */
  bhs.core_orbit_alpha[type][1] = 7.19;  /* alpha_core_2 */
  bhs.core_orbit_c[type][0] = 1.4224;    /* c_core_1     */
  bhs.core_orbit_c[type][1] = -0.4224;   /* c_core_2     */

  /* orbital part */
  eru = 0;
  bhs.alpha[type][eru][1] = 11.13;
  bhs.alpha[type][eru][2] = 13.29;
  bhs.alpha[type][eru][3] = 16.72;
  bhs.c[type][eru][0] = -3.0282;
  bhs.c[type][eru][1] = 0.5619;
  bhs.c[type][eru][2] = 0.0579;
  bhs.c[type][eru][3] = 0.1023;
  bhs.c[type][eru][4] = 0.0040;
  bhs.c[type][eru][5] = -0.0187;

  eru = 1;
  bhs.alpha[type][eru][1] = 9.31;
  bhs.alpha[type][eru][2] = 10.24;
  bhs.alpha[type][eru][3] = 26.07;
  bhs.c[type][eru][0] = 0.3311;
  bhs.c[type][eru][1] = 0.1360;
  bhs.c[type][eru][2] = -0.1867;
  bhs.c[type][eru][3] = 0.1150;
  bhs.c[type][eru][4] = -0.0504;
  bhs.c[type][eru][5] = 0.0051;

  eru = 2;
  bhs.alpha[type][eru][1] = 5.87;
  bhs.alpha[type][eru][2] = 7.12;
  bhs.alpha[type][eru][3] = 8.05;
  bhs.c[type][eru][0] = -1.2035;
  bhs.c[type][eru][1] = -0.7542;
  bhs.c[type][eru][2] = -0.1158;
  bhs.c[type][eru][3] = -0.0806;
  bhs.c[type][eru][4] =  0.0128;
  bhs.c[type][eru][5] =  0.0274;
  /*---------------------------------------------------------------------*/

  for(type=0;type<atom.ntypes;type++) {
    init_bhs_c2a(type); /* evaluate the Ai cofficients */
  }

}

