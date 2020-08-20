#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

/*************************************************************/
/*  [ setting note ]                                         */
/*  Si diamond with H                                        */
/*  one unit cell with 8 Si atoms and 1 H atom               */
/*                                 1999. Dec. project IEMD   */
/*************************************************************/

void  set_param(void)
{
  wv.E_cut  = 5.0;   /* plane wave cut off energy [Hartree]               */
  wv.E_cuti = 2.0;   /* cutoff energy at initial stage of diagonalization */
                     /* E_cut > E_cuti                                    */

  /*--- approach to the Born-Oppenheimer potential surface ---------------*/
  cpq.dt = 4.0;            /* = delta t [a.u.]                */
  cpq.myu = 200.0;	   /* pseudo mass of electron [a.u.]  */
  cpq.residual = 1.0e-5;   /* allowed residual                */
  cpq.max_qstep = 10000;   /* MAX quench step                 */
  cpq.quench = 0.5;	   /* quench rate for wave functios   */
  cpq.qstep_interval = 1;  /* interval for the quenching      */
  cpq.dt2myu = cpq.dt*cpq.dt/cpq.myu;

  /*--- used at CP MD loop -----------------------------------------------*/
  cpa.dt = 4.1336;         /* = delta t [a.u.]                */
  cpa.myu = 200.0;	   /* pseudo mass [a.u.]              */
  cpa.residual = 1.0e-4;   /* allowed residual                */
  cpa.max_step = 30000;	   /* MAX MD step                     */
  cpa.quench = 1.0;	   /* quench rate for atoms motion    */
  cpa.qstep_interval = 1;  /* interval for the quenching      */
  cpa.dt2myu = cpa.dt*cpa.dt/cpa.myu;
  /*----------------------------------------------------------------------*/

  /*--- control temperature ----------------------------------------------*/
  ctl.temp = 1200.0;
  ctl.step = 10;
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

  atom.ntypes = 2;  /* total types of atoms = 1 */

  type = 0;  /*--- [Si] ----------------------------------------------------*/
  atom.m[type]        = 28.0855/(AVOGAD*ME);  /* atomic mass of Si [au]     */
  atom.atomic_n[type] = 14;                   /* atomic number Si = 14      */
  atom.Zv[type]       = 4.0;                  /* valence charge Zv = 4      */
  atom.nvband[type]   = (int)atom.Zv[type];   /* number of valence band     */
                             /*-- check : atom.nvband != atom.Zv in general */
  atom.norbit[type]   = 3;                    /* number of AO : L = 0, 1, 2 */
  atom.ns_orbit[type] = 0;                    /* spin-orbital = none [Si]   */
  /*------------------------------------------------------------------------*/

  type = 1;  /*--- [H] -----------------------------------------------------*/
  atom.m[type]        = 1.00794/(AVOGAD*ME);  /* atomic mass of H [au]      */
  atom.atomic_n[type] = 1;                    /* atomic number H = 1        */
  atom.Zv[type]       = 1.0;                  /* valence charge Zv = 1      */
  atom.nvband[type]   = (int)atom.Zv[type];   /* number of valence band     */
                             /*-- check : atom.nvband != atom.Zv in general */
  atom.norbit[type]   = 3;                    /* number of AO : L = 0, 1, 2 */
  atom.ns_orbit[type] = 0;                    /* spin-orbital = none [H]    */
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
  short	i=0, x, y, z;

  cell.Ax = 5.43*au.r;  /* lattice constant [Si diamond]     */
  cell.Ay = 5.43*au.r;  /* [angstrom] -> [atomic unit]       */          
  cell.Az = 5.43*au.r;
  cell.nx = 1;          /* number of unit cells in each direction */
  cell.ny = 1;
  cell.nz = 1;

  cell.Lx = cell.Ax * (double)cell.nx;   /* basic MD box size */
  cell.Ly = cell.Ay * (double)cell.ny;
  cell.Lz = cell.Az * (double)cell.nz;

  cell.vol = cell.Lx*cell.Ly*cell.Lz;    /* volume of the MD box */ 

  atom.N = 8*cell.nx*cell.ny*cell.nz;    /* total number of atoms [Si] */

  /* normalized coordinates of atoms */
  /*
  cell.rx[0] = 0.00;     cell.rx[4] = 0.50;
  cell.ry[0] = 0.00;     cell.ry[4] = 0.00;
  cell.rz[0] = 0.00;     cell.rz[4] = 0.50;
		                           
  cell.rx[1] = 0.25;     cell.rx[5] = 0.75; 
  cell.ry[1] = 0.25;     cell.ry[5] = 0.75; 
  cell.rz[1] = 0.25;     cell.rz[5] = 0.25; 
		                           
  cell.rx[2] = 0.50;     cell.rx[6] = 0.75;
  cell.ry[2] = 0.50;     cell.ry[6] = 0.25;
  cell.rz[2] = 0.00;     cell.rz[6] = 0.75;
  		                           
  cell.rx[3] = 0.00;     cell.rx[7] = 0.25;
  cell.ry[3] = 0.50;     cell.ry[7] = 0.75;
  cell.rz[3] = 0.50;     cell.rz[7] = 0.75;
  */

  for(z=0; z<cell.nz*2; z++)
    for(y=0; y<cell.ny*2; y++)
      for(x=0; x<cell.nx*2; x++){
	if((x+y+z)%2==0){
	  cell.rx[i] = ((double)x)/2.0;
	  cell.ry[i] = ((double)y)/2.0;
	  cell.rz[i] = ((double)z)/2.0;
	  i++;
	  cell.rx[i] = cell.rx[i-1] + 0.25;
	  cell.ry[i] = cell.ry[i-1] + 0.25;
	  cell.rz[i] = cell.rz[i-1] + 0.25;
	  i++;
	}
      }

  for(i=0;i<atom.N;i++){
    atom.type[i] = 0;
  }

  /*------ adding the [H] atom !!! ------*/
  atom.N += 1; 
  cell.rx[atom.N-1] = 0.10;
  cell.ry[atom.N-1] = 0.10;
  cell.rz[atom.N-1] = 0.10;
  atom.type[atom.N-1] = 1;
  /*-------------------------------------*/

  for(i=0; i<atom.N; i++) {
    cell.rx[i] *= cell.Ax;
    cell.ry[i] *= cell.Ay;
    cell.rz[i] *= cell.Az;
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

  /*--- Si BHS pseudo potentail parameters ------------------------------*/

  type = 0; /* Si */

  /* core part */
  bhs.core_orbit_alpha[type][0] = 2.16;  /* alpha_core_1 */
  bhs.core_orbit_alpha[type][1] = 0.86;  /* alpha_core_2 */
  bhs.core_orbit_c[type][0] = 1.6054;    /* c_core_1     */
  bhs.core_orbit_c[type][1] = -0.6054;   /* c_core_2     */

  /* orbital part */
  eru = 0;
  bhs.alpha[type][eru][1] = 2.48;
  bhs.alpha[type][eru][2] = 2.81;
  bhs.alpha[type][eru][3] = 3.09;
  bhs.c[type][eru][0] = -3.0575;  bhs.c[type][eru][1] =  0.8096;
  bhs.c[type][eru][2] =  0.0012;  bhs.c[type][eru][3] =  0.0511;
  bhs.c[type][eru][4] = -0.0217;  bhs.c[type][eru][5] = -0.0128;

  eru = 1;
  bhs.alpha[type][eru][1] = 1.24;
  bhs.alpha[type][eru][2] = 1.60;
  bhs.alpha[type][eru][3] = 2.12;
  bhs.c[type][eru][0] = -1.7966;  bhs.c[type][eru][1] = -0.0986;
  bhs.c[type][eru][2] =  0.0424;  bhs.c[type][eru][3] =  0.0284;
  bhs.c[type][eru][4] = -0.0030;  bhs.c[type][eru][5] = -0.0039;

  eru = 2;
  bhs.alpha[type][eru][1] = 1.89;
  bhs.alpha[type][eru][2] = 2.22;
  bhs.alpha[type][eru][3] = 2.48;
  bhs.c[type][eru][0] = -0.1817;  bhs.c[type][eru][1] = -0.5634;
  bhs.c[type][eru][2] = -0.0944;  bhs.c[type][eru][3] = -0.2168;
  bhs.c[type][eru][4] =  0.0215;  bhs.c[type][eru][5] =  0.0588;
  /*---------------------------------------------------------------------*/

  /*--- H BHS pseudo potentail parameters -------------------------------*/

  type = 1; /* H */

  /* core part */
  bhs.core_orbit_alpha[type][0] = 16.22;  /* alpha_core_1 */
  bhs.core_orbit_alpha[type][1] = 5.55;   /* alpha_core_2 */
  bhs.core_orbit_c[type][0] = 1.1924;     /* c_core_1     */
  bhs.core_orbit_c[type][1] = -0.1924;    /* c_core_2     */

  /* orbital part */
  eru = 0;
  bhs.alpha[type][eru][1] = 17.08;
  bhs.alpha[type][eru][2] = 23.54;
  bhs.alpha[type][eru][3] = 25.42;
  bhs.c[type][eru][0] =  0.0950;  bhs.c[type][eru][1] = -0.0842;
  bhs.c[type][eru][2] = -0.0443;  bhs.c[type][eru][3] = -0.0519;
  bhs.c[type][eru][4] =  0.0084;  bhs.c[type][eru][5] =  0.0122;

  eru = 1;
  bhs.alpha[type][eru][1] = 1.71;
  bhs.alpha[type][eru][2] = 3.86;
  bhs.alpha[type][eru][3] = 8.08;
  bhs.c[type][eru][0] = -0.2475;  bhs.c[type][eru][1] = -0.2727;
  bhs.c[type][eru][2] = -0.0301;  bhs.c[type][eru][3] = -0.0094;
  bhs.c[type][eru][4] = -0.0291;  bhs.c[type][eru][5] =  0.0018;

  eru = 2;
  bhs.alpha[type][eru][1] = 2.44;
  bhs.alpha[type][eru][2] = 3.30;
  bhs.alpha[type][eru][3] = 4.47;
  bhs.c[type][eru][0] = -0.6887;  bhs.c[type][eru][1] =  0.0913;
  bhs.c[type][eru][2] = -0.1809;  bhs.c[type][eru][3] = -0.0881;
  bhs.c[type][eru][4] = -0.0541;  bhs.c[type][eru][5] =  0.0058;
  /*---------------------------------------------------------------------*/

  for(type=0;type<atom.ntypes;type++) {
    init_bhs_c2a(type); /* evaluate the Ai cofficients */
  }

}

