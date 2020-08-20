#include  <stdio.h>
#include  <sys/types.h>
#include  <time.h>

#include  "prototypes.h"
#include  "physical_constants.h"
#include  "function_macros.h"
#include  "size.h"

/*------------------ global variables -------------------------------------*/
typedef	struct {
  int	  nx;    /* number of unit cells in x-direction          */
  int	  ny;    /*                         y-direction          */
  int	  nz;    /*                         z-direction          */
  double  rx[MAX_ATOMS]; /* positions of atoms in x-coordination */
  double  ry[MAX_ATOMS]; /*                       y-coordination */
  double  rz[MAX_ATOMS]; /*                       z-coordination */
  double  Ax;	 /* initial lattice constant of x-direction      */
  double  Ay;	 /*                             y-direction      */
  double  Az;	 /*                             z-direction      */
  double  Lx;	 /* basic MD cell size = sys.Ax * (double)sys.nx */
  double  Ly;	 /*                      sys.Ay * (double)sys.ny */
  double  Lz;	 /*                      sys.Az * (double)sys.nz */
  double  vol;	 /* volume of the MD cell [a.u.^3]               */
  double  vx[MAX_ATOMS];  /* velocity of atoms in x-coordination */
  double  vy[MAX_ATOMS];  /*                      y-coordination */
  double  vz[MAX_ATOMS];  /*                      z-coordination */
} MD_cell;

typedef struct{  /* BHS pseudo potentail */
  double core_orbit_alpha[MAX_ATOM_TYPE][2];
  double core_orbit_c[MAX_ATOM_TYPE][2];
  double alpha[MAX_ATOM_TYPE][4][4];
  double c[MAX_ATOM_TYPE][4][6];
  double A[MAX_ATOM_TYPE][4][7];

  double wi[10];       /* = weighting for the Gauss-Hermite integral     */
  double xi[10];       /* = nodes for the Gauss-Hermite integral         */
  double c_y[10][10];  /* = coefficients for the spherical harmonics     */

  double bsl[4][MAX_ATOM_TYPE][BESSEL_MAX][BESSEL_MAX];
  int    id[MAX_WAVE_VECTOR2];
} BHS_pseudo_potential_set;

typedef	struct {
  int	               N; /* total number of atoms in the system */
  int	 type[MAX_ATOMS]; /* type of atoms : 0, 1, 2,...         */
  double    m[MAX_ATOMS]; /* atomic mass for each atom types     */
  double   Zv[MAX_ATOMS]; /* valence charge for each atom types  */

  int    atomic_n[MAX_ATOMS]; /* atomic number                            */
  int                 ntypes; /* total number of atom types in the system */

  int    norbit[MAX_ATOMS];   /* number of AO                  */
  int    ns_orbit[MAX_ATOMS]; /* number of spin-AO             */
  int    nvband[MAX_ATOMS];   /* = (int)Zv : number of v-bands */
  int    nvband_all;          /* total number of v-bands       */
  double fermi[MAX_BANDS];    /* fermi distribution            */
} atom_property;

typedef struct{    /* plane wave basis set */
  double E_cut;    /* cutoff energy for the plane wave                    */
  double E_cuti;   /* cutoff energy at initial stage of diagonalization   */
                   /* E_cut > E_cutini                                    */

  int ngx[MAX_WAVE_VECTOR2];     /* reciprocal lattice points             */
  int ngy[MAX_WAVE_VECTOR2];
  int ngz[MAX_WAVE_VECTOR2];
  double kgx[MAX_WAVE_VECTOR2]; 
  double kgy[MAX_WAVE_VECTOR2]; 
  double kgz[MAX_WAVE_VECTOR2];
  double norm2[MAX_WAVE_VECTOR2];              /* index for the FFT3D    */
  int    index[MAX_WAVE_VECTOR2];              /* index for the FFT3D    */
  double fftdat[2*FFTPTSX*FFTPTSY*FFTPTSZ];    /* data (real, imaginary) */

  /* plane wave expansion coefficients [real] [imaginary] */
  double cg1r[MAX_BANDS][MAX_WAVE_VECTOR2];  /* C(t-dt) */
  double cg1i[MAX_BANDS][MAX_WAVE_VECTOR2];

  double cg2r[MAX_BANDS][MAX_WAVE_VECTOR2];  /* C(t)    */
  double cg2i[MAX_BANDS][MAX_WAVE_VECTOR2];

  double cg3r[MAX_BANDS][MAX_WAVE_VECTOR2];  /* C(t+dt) */
  double cg3i[MAX_BANDS][MAX_WAVE_VECTOR2];

  double hcgr[MAX_BANDS][MAX_WAVE_VECTOR2];  /* Hps*C(t) */
  double hcgi[MAX_BANDS][MAX_WAVE_VECTOR2];

  double hvxr[MAX_BANDS][MAX_WAVE_VECTOR2];  /* Hxc*C(t) */
  double hvxi[MAX_BANDS][MAX_WAVE_VECTOR2];

  double rhcr[MAX_BANDS][MAX_WAVE_VECTOR2];  /* H*C(t) */
  double rhci[MAX_BANDS][MAX_WAVE_VECTOR2];  /* H*C(t) */

  int    nplw;     /* total number of reciprocal lattice vectors  */
  int    nplwin;   /*                          for initial stage  */ 

  /* Ryckaert method (next.c) */
  double ncgr[MAX_BANDS][MAX_WAVE_VECTOR2];  /* non orthogonal C(t+dt) */
  double ncgi[MAX_BANDS][MAX_WAVE_VECTOR2];

} plane_wave_set;

typedef struct{
  double dt;       /* = delta t [a.u.]    */
  double myu;      /* pseudo mass [a.u.]  */
  double dt2myu;   /* = dt^2/myu          */
  double residual;
  int    step;
  int    max_qstep;
  int    qstep_interval;
  double quench;
} Car_Parrinello_Quench;

typedef struct{
  double dt;       /* = delta t [a.u.]    */
  double myu;      /* pseudo mass [a.u.]  */
  double dt2myu;   /* = dt^2/myu          */
  double residual;
  int    step;
  int    max_step;
  int    qstep_interval;
  double quench;
} Car_Parrinello_Atom;

typedef struct{   /* Hamiltonian matrix */
  double mat_re[MAX_WAVE_VECTOR1*(MAX_WAVE_VECTOR1+1)/2];
  double mat_im[MAX_WAVE_VECTOR1*(MAX_WAVE_VECTOR1+1)/2];
  double egv[MAX_BANDS];  /* store the engenvalues */
} hamiltonian_set;

typedef struct{   /* atomic unit converters */
  double hbar;
  double mH;
  double mRy;
  double e2H;
  double e2Ry;
  double r;

  double timeH;
  double timeRy;
  double EH2J;
  double ERy2J;
} atomic_unit;

typedef struct{   /* eigenvalues and eigenvectors */
  double d[MAX_WAVE_VECTOR1];
  double d1[MAX_WAVE_VECTOR1];
  double d2[MAX_WAVE_VECTOR1];
  double d3[MAX_WAVE_VECTOR1];
  double a[3][MAX_WAVE_VECTOR1];
  double wk1[MAX_WAVE_VECTOR1];
  double wk2[5][MAX_WAVE_VECTOR1];
  double w[MAX_WAVE_VECTOR1];
  int    iwk[MAX_WAVE_VECTOR1];
  int    iflg[MAX_BANDS];
} eigen_set; 

typedef struct{   /* working valiables */
  double Vps;
  double Vpsmat[MAX_ATOM_TYPE][MAX_WAVE_VECTOR2][MAX_WAVE_VECTOR2];
  double fccx[MAX_ATOMS]; /* core-core EWALD force */
  double fccy[MAX_ATOMS];
  double fccz[MAX_ATOMS];

  double fecx[MAX_ATOMS]; /* electron-core force (Hellmann-Feynman) */
  double fecy[MAX_ATOMS];
  double fecz[MAX_ATOMS];

  double fx[MAX_ATOMS]; /* force = EWALD + HF */
  double fy[MAX_ATOMS];
  double fz[MAX_ATOMS];

  double fx0[MAX_ATOMS]; /* old force */
  double fy0[MAX_ATOMS];
  double fz0[MAX_ATOMS];
} work_set;

typedef struct{  /* for FFT */
  int nx;
  int ny;
  int nz;
  int mesh;
  double rho[FFTPTSX*FFTPTSY*FFTPTSZ]; /* charge density */
} fft_set;

typedef struct{  /* for EWALD */
  double alpha;  /* alpha parameter setting */
  double alpha2;
  double hmax2;  /* = hmax^2 : cut off */
  double kx[EWALD_VECTORS]; /* reciprocal lattice vectors */
  double ky[EWALD_VECTORS];
  double kz[EWALD_VECTORS];
  double knorm2[EWALD_VECTORS];
  double pot;
  int    N;
} ewald_calc_set;

typedef struct{
  double V[FFTPTSX*FFTPTSY*FFTPTSZ];
  double E[FFTPTSX*FFTPTSY*FFTPTSZ];
  double r[FFTPTSX*FFTPTSY*FFTPTSZ];    /* real part */
  double i[FFTPTSX*FFTPTSY*FFTPTSZ];    /* imaginary part */
  double Vhxc[FFTPTSX*FFTPTSY*FFTPTSZ]; /* Hartree + Vxc  potential */
} exchange_correlation_set;

typedef struct{
  double kin, ce, ee, xcen, el;
  double ckin, cct;
  double totf, totr;
  double fker;
} total_energy_set;

typedef struct{
  double egvl[MAX_BANDS]; /* approximated eigen value */
  int    oder[MAX_BANDS];

  double Lr[MAX_BANDS];   /* Legendre Multiplier */
  double Li[MAX_BANDS];

  /* Ryckaert method (next.c) */
  double X0r[MAX_BANDS][MAX_BANDS];
  double X0i[MAX_BANDS][MAX_BANDS];

  double residual[MAX_BANDS];
  double residual_max;
} Car_Parrinello_Next;

typedef struct{
  int step;
  double temp;
  double H2K;
} control_temperature;

/*------------------- declaration for the structures ----------------------*/
  BHS_pseudo_potential_set bhs;
  plane_wave_set           wv;
  hamiltonian_set          hm;
  atomic_unit              au;
  atom_property            atom;
  MD_cell                  cell;
  eigen_set                eg;
  work_set                 wk;
  fft_set                  ft;
  ewald_calc_set           ew;
  exchange_correlation_set hxc;
  total_energy_set         E;
  Car_Parrinello_Quench    cpq;
  Car_Parrinello_Atom      cpa;
  Car_Parrinello_Next      cp;
  control_temperature      ctl;
/*-------------------------------------------------------------------------*/
