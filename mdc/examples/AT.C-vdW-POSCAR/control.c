#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*********************************************************************/
/*  Carbon nanotube (21,21) with Abell Tersoff potential             */
/*********************************************************************/
/*****
*  purpose: system initialisation
*  all valiables are specified in a MD calculational unit.
*  see also [manual/units.IEMD]
*****/
FILE *fp_poscar;

void read_poscar(void)
{
  int i=0;

  if((fp_poscar = fopen("./POSCAR","r"))==NULL){
    printf("cannot open poscar. Abort\n");
    exit(EXIT_FAILURE);
  }
  fgets(vsp.char_dumy, sizeof(vsp.char_dumy), fp_poscar); /* read dumy char */
  fscanf(fp_poscar, "%lf", &vsp.eru);

  /* super cell definition */
  fscanf(fp_poscar, "%lf", &vsp.eru_xx); vsp.eru_xx*=vsp.eru;  
  fscanf(fp_poscar, "%lf", &vsp.eru_xy); vsp.eru_xy*=vsp.eru;
  fscanf(fp_poscar, "%lf", &vsp.eru_xz); vsp.eru_xz*=vsp.eru;

  fscanf(fp_poscar, "%lf", &vsp.eru_yx); vsp.eru_yx*=vsp.eru;
  fscanf(fp_poscar, "%lf", &vsp.eru_yy); vsp.eru_yy*=vsp.eru;
  fscanf(fp_poscar, "%lf", &vsp.eru_yz); vsp.eru_yz*=vsp.eru;

  fscanf(fp_poscar, "%lf", &vsp.eru_zx); vsp.eru_zx*=vsp.eru;
  fscanf(fp_poscar, "%lf", &vsp.eru_zy); vsp.eru_zy*=vsp.eru;
  fscanf(fp_poscar, "%lf", &vsp.eru_zz); vsp.eru_zz*=vsp.eru;

  fscanf(fp_poscar, "%d", &vsp.n_atom);  /* total number of atoms */
  printf("%d", vsp.n_atom);

  fgets(vsp.char_dumy, sizeof(vsp.char_dumy), fp_poscar); /* read dumy char */
  fgets(vsp.char_dumy, sizeof(vsp.char_dumy), fp_poscar); /* read dumy char */

  /* particle positions */	
  vsp.rx = (double *)calloc(vsp.n_atom, sizeof(double)); 
  vsp.ry = (double *)calloc(vsp.n_atom, sizeof(double)); 
  vsp.rz = (double *)calloc(vsp.n_atom, sizeof(double)); 


  for(i=0;i<vsp.n_atom;i++){  /* only for the orthogonal cell */
    fscanf(fp_poscar, "%lf", &vsp.rx[i]); vsp.rx[i]*=vsp.eru_xx;
    fscanf(fp_poscar, "%lf", &vsp.ry[i]); vsp.ry[i]*=vsp.eru_yy;
    fscanf(fp_poscar, "%lf", &vsp.rz[i]); vsp.rz[i]*=vsp.eru_zz;
    /*    printf(" %lf %lf %lf\n",vsp.rx[i],vsp.ry[i],vsp.rz[i]); */
  }
  fclose(fp_poscar);
}

void get_control_param(void)
{

  read_poscar(); /* read the POSCAR */
                 /* POSCAR defines the atomic positions in the system.   */
                 /* The file is a input file for the ab-initio code VASP */

  ctl.calc_max = 50000;     /* maximum MD time step                     */
  ctl.delta_time_fs = 0.3;  /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 20.0;            /* Temperature setting [K]                  */
  ctl.t_control_step = 100;   /* scale at every ctl.t_control_step steps  */
  ctl.set_press_GPa_X = 0.01; /* Pressure setting  [GPa]                  */
  ctl.set_press_GPa_Y = 0.01; 
  ctl.set_press_GPa_Z = 0.01; 
  ctl.p_control_step  =  20;
  sys.pres_scalling_factor = 0.5;

  sys.Ax = vsp.eru_xx;   /* super_cell */
  sys.Ay = vsp.eru_yy;
  sys.Az = vsp.eru_zz;
  sys.nx = 1;	   /* number of unit cells in x-direction                 */
  sys.ny = 1;      /* note: currently only nx=ny=nz=1 of C60 is supported */
  sys.nz = 1;

  ctl.natoms_in_unit_cell = vsp.n_atom;  /* number of atoms */
  ctl.natoms_in_mol_unit  = vsp.n_atom; 
  ctl.kinds_of_ions       = 1;     /* C */

  /* alpha setting for EWALD calculation */
  sys.a1 = 0.0;     /* to SKIP the EWALD calculation */
  sys.hm = 20;      /* has no meaning for AT-C */
  sys.radius = 3.0; /* [A] AT potential set */  

}

/*****
  identify_ion:
  every atom is "C" in C crystal with diamond structure 
*****/
void	identify_ion(void)

{
  short	x, y, z, i=0;

  for(z=0; z<sys.nz*2; z++)
    for(y=0; y<sys.ny*2; y++)
      for(x=0; x<sys.nx*2; x++)
	sys.ion[i++] = 0;
}


/*****
  void set_loc(void)
*****/
void	set_loc(void)
{
  int i;

  for(i=0;i<vsp.n_atom;i++){  /* only for the orthogonal cell */
    sys.rx[i] = vsp.rx[i];
    sys.ry[i] = vsp.ry[i];
    sys.rz[i] = vsp.rz[i];
  }


}

void set_potential(void)
{
  double m_C;

  m_C = 12.011 /(6.0221367 * 1.6605402 / 10.0); /* [m_u] */

  ion.m[0] = m_C;

  /* for C system : taken from PRB 76, 195447 (2007) */
  at.A = 1.2067090e3 * 1.60217733e-1;  /* [eV -> en] */
  at.B = 3.1596460e2 * 1.60217733e-1;  /* [eV -> en] */

  /* ------------ A.T Potential Set --------------*/
  at.A = 1.3936e3 * 1.60217733e-1;  /* [eV -> en] */
  at.B = 3.4670e2 * 1.60217733e-1;  /* [eV -> en] */

  at.lam = 3.4879;
  at.mu = 2.2119;
  at.beta = 1.5724e-7;
  at.n = 7.2751e-1;
  at.c = 3.8049e4;
  at.d = 4.384;
  at.h = -0.57058;
  at.R = 1.8;
  at.S = 2.1;               /* = sys.radius = cut-off radius */
  at.x = 1.0;

  /* vdW force for the graphite : PRB 68, 035425 (2003) */
  /* U(r)=4 eps[(sig/r)^12 - (sig/r)^6]                 */

  vdW.eps        = 0.003* 1.60217733e-1;  /* [eV -> en] */
  vdW.sig        = 3.4;                   /* [A]        */
  vdW.sig4       = vdW.sig*vdW.sig*vdW.sig*vdW.sig;
  vdW.sig8       = vdW.sig4*vdW.sig4;
  vdW.sig14      = vdW.sig8*vdW.sig4*vdW.sig*vdW.sig;
  vdW.eps_sig2   = 48.0*vdW.eps/(vdW.sig*vdW.sig);

  vdW.r1    = 3.4; /* vdW force will work only in the region from r1 to r2 */
  vdW.r2    = 4.5;  /* note that the value r2 must be less than the system  */
                    /* size to aboid the undesirable boundary effect        */
}

void mk_table(void)
{
  /*  Carbon (diamond structure) with Abell Tersoff potential      */
  /*  this function do nothing.  see "real.c".                     */
}






