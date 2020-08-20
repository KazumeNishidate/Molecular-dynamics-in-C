/*********************************************************************/
/*  Carbon (diamond structure) with J.Tersoff potential              */
/*********************************************************************/

typedef struct{            /* to read POSCAR file */
  double  *rx, *ry, *rz;   /* positions           */
  double  x, y, z;

  char char_dumy[100];
  int  n_atom;
  double eru;
  double eru_xx;
  double eru_xy;
  double eru_xz;
  double eru_yx;
  double eru_yy;
  double eru_yz;
  double eru_zx;
  double eru_zy;
  double eru_zz;
} vasp;

typedef struct{
  double A;
  double B;
  double lam;
  double mu;
  double beta;
  double n;
  double c;
  double d;
  double h;
  double R;
  double S;
  double x;

} potential_at_set;

  /* vdW force for the graphite : PRB 68, 035425 (2003) */
  /* U(r)=4 eps[(sig/r)^12 - (sig/r)^6] */
typedef struct{ 
  double eps;   /* 0.003 eV */
  double sig;   /* 3.4 A    */  
  double sig4;
  double sig8;
  double sig14;
  double eps_sig2; 
  double r1;    /* cutoff 1 */ 
  double r2;    /* cutoff 2 */
} potential_vdW;

/*------------------- declaration for the structures ----------------------*/
  potential_at_set at;
  vasp vsp;
  potential_vdW  vdW;

