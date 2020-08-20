#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "md.h"
#include  "prototypes.h"

/*********************************************************************/
/* TBMD for silicon                                                  */
/* [Ref.]                                                            */
/*        I.Kwon, R.Biswas, C.Z.Wang, K.M.Ho, and C.M.Soukoulis      */
/*        ``Transferable tight-binding models for silicon'',         */
/*        Physical Review B, Vol.49, Num.11, 7242 (1994).            */
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
  delta_time_fs = 1.07;  /* 1 [fs] = 0.001 [ps] = 1.0 X 10^(-15) [s] */

  ctl.temp = 300.0;       /* Temperature setting [K] */
  ctl.t_control_step = 5;

  set_press_GPa_X = 0.01;  /* Pressure setting  [GPa] */
  set_press_GPa_Y = 0.01; 
  set_press_GPa_Z = 0.01; 
  ctl.p_control_step = 20;
  sys.pres_scalling_factor = 0.5;

  /* check */
  sys.Ax = 5.41959;     /* initial lattice constant */
  sys.Ay = 5.41959;     /* initial lattice constant */
  sys.Az = 5.41959;     /* initial lattice constant */

  sys.nx = 1;	   /* number of unit cells in x-direction */
  sys.ny = 1;
  sys.nz = 1;

  /* check */
  /*--- Diamond structure ---------------------------------*/
  sys.N = 8*sys.nx*sys.ny*sys.nz;  /* total number of ions */
  cluster.type = 0;            
  /*-------------------------------------------------------*/  

  /*--- Silicon cluster ---------------------------------------------*/
  /* sys.N | cluster.type |                    | cohesive energy[eV] */
  /*=================================================================*/
  /*   2   |      1       | dimer              |        1.60         */
  /*-----------------------------------------------------------------*/
  /*   3   |      1       | isosceles-triangle |        2.51         */
  /*       |      2       | linear-chain       |        2.21         */
  /*-----------------------------------------------------------------*/
  /*   4   |      1       | rhombus            |        3.21         */
  /*       |      2       | tetrahedron        |        2.83         */
  /*-----------------------------------------------------------------*/
  /*   5   |      1       | trigonal-bipyramid |        3.18         */
  /*       |      2       | square-pyramid     |        3.13         */
  /*       |      3       | pentagon ?         |        3.04         */
  /*-----------------------------------------------------------------*/
  /*-------------------------*/  
  /*
  sys.nx = sys.ny = sys.nz;
  sys.N = 2;
  cluster.type = 1;
  */
  /*-------------------------*/  

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
  /* TBMD specific part [Si]                                         */
  /*=================================================================*/
  /*-------------  H_tb ---------------------------------------------*/

  /* [eV] -> [Energy] unit converter                          */
  /* 1 [eV] = 1.60217733*10^(-19) [J] -> 0.160217733 [energy] */
  htb.eV2E = 0.160217733;

  /* energy parameters for Htb [eV] -> Htb [Energy] */
  htb.e_0 = 8.7393204 *htb.eV2E;  /* [Energy] */
  htb.e_s = -5.25 *htb.eV2E;  /* [Energy] */
  htb.e_p =  1.20 *htb.eV2E;  /* [Energy] */
  htb.v_A = -2.038 *htb.eV2E;  /* [Energy] */      /* ss sigma --> A */
  htb.v_B =  1.745 *htb.eV2E;  /* [Energy] */      /* sp sigma --> B */
  htb.v_C =  2.75  *htb.eV2E;  /* [Energy] */      /* pp sigma --> C */
  htb.v_D = -1.075 *htb.eV2E;  /* [Energy] */      /* pp pi    --> D */

  /* parameters for function s(r) [1] */
  htb.s_n = 2.0;  /* [dimensionless] */
  htb.s_nc_A = 9.5;  /* [dimensionless] */
  htb.s_nc_B = 8.5;  /* [dimensionless] */
  htb.s_nc_C = 7.5;  /* [dimensionless] */
  htb.s_nc_D = 7.5;  /* [dimensionless] */
  htb.s_rc_A = 3.4;   /* [A] */
  htb.s_rc_B = 3.55;  /* [A] */
  htb.s_rc_C = 3.7;   /* [A] */
  htb.s_rc_D = 3.7;   /* [A] */
  htb.s_r0   = 2.360352;  /* [A] */
  htb.s_r1   = 4.0;   /* cutoff-1   [A]   */ 
  htb.s_r1_2 = 16.0;  /* cutoff-1^2 [A^2] */ 
  htb.s_rm   = 4.16;     /* cutoff-2   [A]   */
  htb.s_rm_2 = 17.3056;  /* cutoff-2^2 [A^2] */

  /* parameters for function s(r) [2] */
  htb.s_c0A =  0.0000317229;  /* [dimensionless] */
  htb.s_c1A = -0.000721508;   /* [A^(-1)]        */
  htb.s_c2A =  0.00530131;    /* [A^(-2)]        */
  htb.s_c3A = -0.0126942;     /* [A^(-3)]        */

  htb.s_c0B =  0.00149092;
  htb.s_c1B = -0.0182203;
  htb.s_c2B =  0.0530362;
  htb.s_c3B =  0.0162594;

  htb.s_c0C =  0.0103039;
  htb.s_c1C = -0.0744898;
  htb.s_c2C = -0.276367;
  htb.s_c3C =  2.12145;

  htb.s_c0D =  0.0103039;
  htb.s_c1D = -0.0744898;
  htb.s_c2D = -0.276367;
  htb.s_c3D =  2.12145;

  /* x^y = exp(y*log(x))  [y = integer, x><0] or [y=float, x>0]  */
  htb.r0_rc_nc_A = exp(htb.s_nc_A*log(htb.s_r0/htb.s_rc_A));
  htb.r0_rc_nc_B = exp(htb.s_nc_B*log(htb.s_r0/htb.s_rc_B));
  htb.r0_rc_nc_C = exp(htb.s_nc_C*log(htb.s_r0/htb.s_rc_C));
  htb.r0_rc_nc_D = exp(htb.s_nc_D*log(htb.s_r0/htb.s_rc_D));

  /*-------------  Repulsion ----------------------------------------*/
  phi.m  =  6.8755;    /* [dimensionless] */
  phi.mc = 13.017;     /* [dimensionless] */
  phi.dc =  3.66995;   /* [A] */
  phi.r0 = htb.s_r0;   /* [A] */
  phi.r1 = htb.s_r1;

  /* x^y = exp(y*log(x))  [y = integer, x><0] or [y=float, x>0]  */
  phi.r0_dc_mc = exp(phi.mc*log(phi.r0/phi.dc));  /* =(r0/dc)^mc */

  /* parameters for repulsive function phi(r) [2] */
  phi.c0p =  1.878970043092705e-11 *htb.eV2E;  /* [Energy]     */
  phi.c1p = -1.322056910835325e-9  *htb.eV2E;  /* [Energy/A]   */
  phi.c2p =  1.43237933661923e-8   *htb.eV2E;  /* [Energy/A^2] */
  phi.c3p = -4.246818966596628e-8  *htb.eV2E;  /* [Energy/A^3] */

  phi.c1f =  2.1604385    *htb.eV2E;  /* [Energy] */
  phi.c2f = -0.1384393    *htb.eV2E;  /* [Energy] */
  phi.c3f =  5.8398423e-3 *htb.eV2E;  /* [Energy] */
  phi.c4f = -8.0263577e-5 *htb.eV2E;  /* [Energy] */

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
  for(i=0; i<sys.N; i++){sys.ion[i] = 0;}
}

/*****
  void set_roc(void)
  set up the particle positions
*****/
void	set_roc(void)
{
  short	x, y, z, i=0;
  short nx, ny, nz;
  double nnx=0.0, nny=0.0, nnz=0.0;
  double nx3=0.0, ny3=0.0, nz3=0.0;
  double shiftx, shifty, shiftz;
  
  double a, b, c, theta;  /* for cluster */

  nx = sys.nx * 2;
  ny = sys.ny * 2;
  nz = sys.nz * 2;

  if(cluster.type!=0){
    for(i=0; i<sys.N; i++) {
      sys.vx[i] = 0.0;
      sys.vy[i] = 0.0;
      sys.vz[i] = 0.0;
    }
  }
  
/**** 
  X display shift 
****/
  shiftx = (sys.Ax)/2.0;
  shifty = (sys.Ay)/2.0;
  shiftz = (sys.Az)/2.0;
  
/****
  fcc ion
****/
  if(cluster.type==0){  /* diamond structure */
    for(z=0; z<nz; z++)
      for(y=0; y<ny; y++)
	for(x=0; x<nx; x++){
	  
	  if((x+y+z)%2==0){
	    nnx = (double)(x)/2.0;
	    nny = (double)(y)/2.0;
	    nnz = (double)(z)/2.0;
	    sys.rx[i] = sys.Ax * nnx +shiftx;
	    sys.ry[i] = sys.Ay * nny +shifty;
	    sys.rz[i] = sys.Az * nnz +shiftz;
	    i++;
	    
	    nx3 = nnx + 0.25;
	    ny3 = nny + 0.25;
	    nz3 = nnz + 0.25;
	    sys.rx[i] = sys.Ax * nx3 +shiftx;
	    sys.ry[i] = sys.Ay * ny3 +shifty;
	    sys.rz[i] = sys.Az * nz3 +shiftz;
	    i++;
	  }
	}
  }
  
/*-----------------------------------------------------------*/  
  if(sys.N==2 && cluster.type==1){  /* dimer */
    a = 2.45;
    sys.rx[0] = sys.Ax/2.0-a/2.0;  sys.rx[1] = sys.Ax/2.0+a/2.0;
    sys.ry[0] = sys.Ay/2.0;        sys.ry[1] = sys.Ay/2.0; 
    sys.rz[0] = sys.Az/2.0;        sys.rz[1] = sys.Az/2.0; 
  }
/*-----------------------------------------------------------*/  
  if(sys.N==3){
    if(cluster.type==1){  /* isosceles-triangle */
      theta = 74.7*(2.0*PI/360.0);
      a = 2.42*cos(theta/2.0);
      b = 2.42*sin(theta/2.0);
      sys.rx[0] = sys.Ax/2.0-a/2.0;  sys.rx[1] = sys.Ax/2.0+a/2.0;  
      sys.ry[0] = sys.Ay/2.0;        sys.ry[1] = sys.Ay/2.0+b;    
      sys.rz[0] = sys.Az/2.0;        sys.rz[1] = sys.Az/2.0;      

      sys.rx[2] = sys.Ax/2.0+a/2.0;
      sys.ry[2] = sys.Ay/2.0-b;  
      sys.rz[2] = sys.Az/2.0;    
    }
    
    if(cluster.type==2){  /* linear-chain */
      a = 2.41;
      sys.rx[0] = sys.Ax/2.0-a;  sys.rx[1] = sys.Ax/2.0;  
      sys.ry[0] = sys.Ay/2.0;    sys.ry[1] = sys.Ay/2.0;  
      sys.rz[0] = sys.Az/2.0;    sys.rz[1] = sys.Az/2.0;  

      sys.rx[2] = sys.Ax/2.0+a;
      sys.ry[2] = sys.Ay/2.0;  
      sys.rz[2] = sys.Az/2.0;  
    }
  }
/*-----------------------------------------------------------*/  
  if(sys.N==4){
    if(cluster.type==1){  /* rhombus */
      c = 2.48;
      b = 2.56/2.0;
      theta = asin(b/c);
      a = c*cos(theta);
      sys.rx[0] = sys.Ax/2.0-a;  sys.rx[1] = sys.Ax/2.0;      
      sys.ry[0] = sys.Ay/2.0;    sys.ry[1] = sys.Ay/2.0+b;      
      sys.rz[0] = sys.Az/2.0;    sys.rz[1] = sys.Az/2.0;  
      
      sys.rx[2] = sys.Ax/2.0;    sys.rx[3] = sys.Ax/2.0+a; 
      sys.ry[2] = sys.Ay/2.0-b;  sys.ry[3] = sys.Ay/2.0;       
      sys.rz[2] = sys.Az/2.0;    sys.rz[3] = sys.Az/2.0;       
    }
    
    if(cluster.type==2){  /* tetrahedron */
      a = 2.59/2.0;
      b = sqrt(3.0/4.0)*a;
      c = sqrt(2.0/3.0)*a;
      sys.rx[0] = sys.Ax/2.0-b;  sys.rx[1] = sys.Ax/2.0-b;
      sys.ry[0] = sys.Ay/2.0-a;  sys.ry[1] = sys.Ay/2.0+a;              
      sys.rz[0] = sys.Az/2.0-c;  sys.rz[1] = sys.Az/2.0-c;
      
      sys.rx[2] = sys.Ax/2.0+b;  sys.rx[3] = sys.Ax/2.0-b+2.0*a/sqrt(12.0); 
      sys.ry[2] = sys.Ay/2.0;    sys.ry[3] = sys.Ay/2.0;
      sys.rz[2] = sys.Az/2.0-c;  sys.rz[3] = sys.Az/2.0+c; 
    }
  }
/*-----------------------------------------------------------*/  
  if(sys.N==5){
    if(cluster.type==1){  /* trigonal-bipyramid */
      a = 3.59/2.0; 
      b = 2.74/2.0; 
      c = a*sqrt(3.0)/3.0;
      sys.rx[0] = sys.Ax/2.0-c;      sys.rx[1] = sys.Ax/2.0-c;
      sys.ry[0] = sys.Ay/2.0-a;      sys.ry[1] = sys.Ay/2.0+a;
      sys.rz[0] = sys.Az/2.0;        sys.rz[1] = sys.Az/2.0;      

      sys.rx[2] = sys.Ax/2.0+2.0*c;  sys.rx[3] = sys.Ax/2.0;  
      sys.ry[2] = sys.Ay/2.0;        sys.ry[3] = sys.Ay/2.0;
      sys.rz[2] = sys.Az/2.0;        sys.rz[3] = sys.Az/2.0+b;      

      sys.rx[4] = sys.Ax/2.0;
      sys.ry[4] = sys.Ay/2.0;
      sys.rz[4] = sys.Az/2.0-b;
    }

    if(cluster.type==2){  /* square-pyramid */
      a = 2.56/2.0; 
      b = a/sqrt(2.0);
      sys.rx[0] = sys.Ax/2.0-a;  sys.rx[1] = sys.Ax/2.0+a;
      sys.ry[0] = sys.Ay/2.0-a;  sys.ry[1] = sys.Ay/2.0-a;
      sys.rz[0] = sys.Az/2.0-b;  sys.rz[1] = sys.Az/2.0-b;      
      
      sys.rx[2] = sys.Ax/2.0-a;  sys.rx[3] = sys.Ax/2.0+a;
      sys.ry[2] = sys.Ay/2.0+a;  sys.ry[3] = sys.Ay/2.0+a;
      sys.rz[2] = sys.Az/2.0-b;  sys.rz[3] = sys.Az/2.0-b;      

      sys.rx[4] = sys.Ax/2.0;
      sys.ry[4] = sys.Ay/2.0;
      sys.rz[4] = sys.Az/2.0+b;
    }
    
    if(cluster.type==3){  /* pentagon */  
    }
  }
/*-----------------------------------------------------------*/
}

void set_potential(void)
{
  ion.m[0] = 28.0855/(6.0221367*1.6605402/10.0);  /* mass of Si atom [m_u] */
}
