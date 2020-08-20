#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "md.h"
#include "potential.h"
#include "prototypes.h"

int    max_common_measure(int, int);
double norm_of_T(int, int);
void   atoms_in_chain(void);
int    calc_atoms( int, int );

int max_common_measure( int n, int m )
{
  if( m == 0 ) return n;
  return max_common_measure( m, (n%m) );
}

double norm_of_T(int n, int m)
{

  int dR, mcm, mcm3;
  double norm;

  if(n < m){
    printf("Your input is strange!!\n");
    printf("n must be larger than m.\n");
    exit(0);
  }

  mcm  = max_common_measure(n, m);
  mcm3 = 3 * mcm;

  if((n - m) % mcm3 == 0){
    dR = mcm3;
  }else{
    dR = mcm;
  }

  norm = 3.0 * tube.bond_length * sqrt((double)(n * n + n * m + m * m));

  return norm / (double)dR;

}

void atoms_in_chain(void)
{

  int i, j, n, num, flag, type;
  double x0, y0, x, y, dx1, dx2, dy1, dy2;
  double phi, psi, theta;
  double si, co, bond_length;
  /*-------*/
  int d;
  double  dx, dy;
  double *temp_point_y, *temp_point_x;
  /*-------*/

  phi = tube.ch_theta - (PI / 6.0);
  psi = tube.ch_theta + (PI / 6.0);

  theta = PI / 2.0 + tube.ch_theta;
  si = sin(theta);
  co = cos(theta);

  dx1 = tube.bond_length * cos(phi);
  dy1 = tube.bond_length * sin(phi);
  dx2 = tube.bond_length * cos(psi);
  dy2 = tube.bond_length * sin(psi);

  /*-----------------------------------------------------------------------*/
  temp_point_y = (double *)calloc((2 * tube.ch_n + tube.ch_m), sizeof(double));
  temp_point_x = (double *)calloc((2 * tube.ch_n + tube.ch_m), sizeof(double));

  n = 1;
  x = 0.0;
  y = tube.T;
  for(i = 0; i < (2 * tube.ch_n + tube.ch_m); i++){
    temp_point_x[i] = x;
    temp_point_y[i] = y;

    if(n == 1){
      x = x + dx1;
      y = y + dy1;
      n = 2;
    }else{
      x = x + dx2;
      y = y + dy2;
      n = 1;
    }
  }
  /*-----------------------------------------------------------------------*/


  y0 = 0.0;
  x0 = 0.0;
  for(i = 0; i < (2 * tube.ch_n + tube.ch_m); i++){

    n = 1; num = 0; flag = 0;
    type = i % 2; y = y0; x = x0;

    while(flag == 0){

      num++;
      if(n == 1){
	bond_length = tube.bond_length;
	if(type == 1) bond_length = 2.0 * tube.bond_length;
	n = 2;
      }else{
	bond_length = 2.0 * tube.bond_length;
        if(type == 1) bond_length = tube.bond_length;
	n = 1;
      }

      y = y + bond_length * si;
      x = x + bond_length * co;
      if(x < 0.0) x += tube.norm_ch;
      if(x > tube.norm_ch) x -= tube.norm_ch;

      for(j = 0; j < (2 * tube.ch_n + tube.ch_m); j++){

	dx = (temp_point_x[j] - x) * 100.0;
	if(dx < 0.0) dx = dx * (-1.0);
	d = (int)dx;
	if(d == 0){

	  dy = (temp_point_y[j] - y) * 100.0;
	  if(dy < 0.0) dy = dy * (-1.0);
	  d = (int)dy;
	  if(d == 0){
	    flag = 1;
	    break;

	  }

	}
      }

    }/* <--- end of while */

    tube.chain[i] = num;

    if(type == 0){
      y0 = y0 + dy1;
      x0 = x0 + dx1;
    }else{
      y0 = y0 + dy2;
      x0 = x0 + dx2;
    }

  }/* <--- end of loop i */

  free(temp_point_x);
  free(temp_point_y);

}

void nanotube_init(int n, int m, double a)
/*************************************************/
/**  (n, m) is a index number of chiral vector  **/
/**                                             **/
/**    a is bond-length between carbon atoms    **/
/*************************************************/
{

  int i;
  double temp = 0.0;

  if( tube.chain != NULL ) tube.chain = NULL;
  if( tube.rx != NULL ) tube.rx = NULL;
  if( tube.ry != NULL ) tube.ry = NULL;
  if( tube.rz != NULL ) tube.rz = NULL;

  tube.chain = (int *)calloc((2 * n + m), sizeof(int));

  tube.ch_n = n;
  tube.ch_m = m;

  tube.norm_ch = 1.732050807 * a * sqrt((double)(n * n + n * m + m * m));

  tube.bond_length = a;

  temp = 1.732050807 * m / ((double)(2 * n + m));
  tube.ch_theta = atan(temp);

  tube.rad = tube.norm_ch / 2.0 / PI;

  tube.T = norm_of_T(n, m);

  atoms_in_chain();

  tube.atoms = 0;
  for(i = 0; i < (2 * n + m); i++){
    tube.atoms += tube.chain[i];
  }

  tube.rx = ( double * )calloc( sys.ny*tube.atoms, sizeof( double ) );
  tube.ry = ( double * )calloc( sys.ny*tube.atoms, sizeof( double ) );
  tube.rz = ( double * )calloc( sys.ny*tube.atoms, sizeof( double ) );

  create_nanotube();
}

void chiral_structure(double x, double y, double z, int j, int type)
 /*********************************************************/
 /**  Point(x, y, z) is a start point of a chiral chain  **/
 /**                                                     **/
 /**       type : 1  --->  o--o----o--o----o......       **/
 /**              2  --->  o----o--o----o--o......       **/
 /*********************************************************/        
{

  static int k = 0;
  int i;
  double temp;
  double phi, theta;
  double d, bond_length;
  double xx, yy, zz;

  xx = x; yy = y; zz = z;

  theta = PI / 2.0 + tube.ch_theta;

  d   = 0.0;
  phi = 0.0;

  for(i = 0; i < tube.chain[j]; i++){
    temp = xx;
    xx   = temp * cos(phi) - zz * sin(phi);
    zz   = temp * sin(phi) + zz * cos(phi);
    yy   = yy + d;

    /** cyclic boundary condition for y-direction **/
    if(yy >= ( tube.T*sys.ny )) yy = yy - sys.Ly;
    if(yy < 0.0)     yy = yy + sys.Ly; 

    tube.rx[k] = xx;
    tube.rz[k] = zz;
    tube.ry[k] = yy;

    k++;

    if((i % 2) != 0){
      bond_length = 2.0 * tube.bond_length;
      if(type == 2) bond_length = tube.bond_length;
    }else{
      bond_length = tube.bond_length;
      if(type == 2) bond_length = 2.0 * tube.bond_length;
    }

    d   = bond_length * sin(theta);
    phi = bond_length * cos(theta) / tube.rad;
  }

}

void armchair(void)
{
  
  int i, j, k = 0, n, sign;
  double temp;
  double d, x, y, z;
  double co, si, phi, psi;
  
  phi = PI / 3.0;
  co = cos(phi);
  si = sin(phi);
  
  psi = tube.bond_length / tube.rad;
  phi = tube.bond_length * co / tube.rad;
  d = tube.bond_length * si;

  n = tube.atoms;

  for(j = 0; j < sys.ny; j++){
    sign = 1;
    x = tube.rad;
    y = 0.0 + tube.T * (double)j;
    z = 0.0;
    for(i = 0; i < n; i++){
      tube.rx[k] = x;
      tube.ry[k] = y;
      tube.rz[k] = z;
    
      if((i % 2) == 0){
	temp = x;
	x = temp * cos(psi) - z * sin(psi);
	z = temp * sin(psi) + z * cos(psi);
	y = y;
      }else{
	temp = x;
	x = temp * cos(phi) - z * sin(phi);
	z = temp * sin(phi) + z * cos(phi);
	y = y + d * sign;
	sign = -sign;
      }

      k++;

    }

  }

}

void create_nanotube( void )
{
  int i, j, n;
  double d, temp;
  double psi, xi;
  double x1, y1, z1;
  double x2, y2, z2;

  if(tube.ch_n == tube.ch_m){
    armchair();
    return;
  }

  for(n = 0; n < sys.ny; n++){

    x1 = tube.rad; z1 = 0.0;
    y1 = (double)n * tube.T;

    i  = 0;
    d  = 0.0;
    xi = 0.0;

    for(j = 0; j < (2 * tube.ch_n + tube.ch_m); j++){

      if((j % 2) == 0){                   /* <--- Type I chain */
        temp = x1;
        x1 = temp * cos(xi) - z1 * sin(xi);
        z1 = temp * sin(xi) + z1 * cos(xi);
        y1 = y1 + d;

        chiral_structure(x1, y1, z1, j, 1);

      }else{                              /* <--- Type II chain */
        psi = PI / 6.0 - tube.ch_theta;
        y2  = y1 - tube.bond_length * sin(psi);
        psi = tube.bond_length * cos(psi) / tube.rad;
        x2  = x1 * cos(psi) - z1 * sin(psi);
        z2  = x1 * sin(psi) + z1 * cos(psi);

        chiral_structure(x2, y2, z2, j, 2);

      }

      if(i == 0){
        i++;
        d  = 1.732050807 * tube.bond_length;
        xi = d * cos(tube.ch_theta) / tube.rad;
        d  = d * sin(tube.ch_theta);
      }

    }/* <--- end of j loop */

  }/* <--- end of n loop */

}

void delete_nanotube( void )
{
    if( tube.chain != NULL ) free( tube.chain );
    if( tube.rx != NULL ) free( tube.rx );
    if( tube.ry != NULL ) free( tube.ry );
    if( tube.rz != NULL ) free( tube.rz );
}

int calc_atoms( int n, int m )
{
    int ret, dR, mcm, mcm3;

    mcm  = max_common_measure(n, m);
    mcm3 = 3 * mcm;

    if((n - m) % mcm3 == 0){
      dR = mcm3;
    }else{
      dR = mcm;
    }
    ret = 4*( n*n + m*m + n*m )/dR; 

    return ret;
}
/*--------------------------------------------------------------------------*/
void multiwall_nanotube( int wall, int n1, int m1, int n2, int m2 )
{
    int i, k = 0, num = 0;
    multi.num_of_wall = wall;

    num += calc_atoms( n1, m1 );
    if( wall == 2 ) num += calc_atoms( n2, m2 );
    multi.atoms = num;
    multi.rx = ( double * )calloc( num*sys.ny, sizeof( double ) );
    multi.ry = ( double * )calloc( num*sys.ny, sizeof( double ) );
    multi.rz = ( double * )calloc( num*sys.ny, sizeof( double ) );

    nanotube_init( n1, m1, 1.42 );
    for( i = 0; i < tube.atoms*sys.ny; i++ ){
        multi.rx[ k ] = tube.rx[ i ];
        multi.ry[ k ] = tube.ry[ i ];
        multi.rz[ k ] = tube.rz[ i ];
        k++;
    }
    if( wall == 2 ){
      nanotube_init( n2, m2, 1.42 );
      for( i = 0; i < tube.atoms*sys.ny; i++ ){
          multi.rx[ k ] = tube.rx[ i ];
          multi.ry[ k ] = tube.ry[ i ];
          multi.rz[ k ] = tube.rz[ i ];
          k++;
      }
    }

    multi.T = tube.T;
    multi.rad = tube.rad;

    delete_nanotube(); 
}

void delete_multi_nanotube( void )
{
    if( multi.rx != NULL ) free( multi.rx );
    if( multi.ry != NULL ) free( multi.ry );
    if( multi.rz != NULL ) free( multi.rz );
}
