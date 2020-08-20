#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*****
  void	look_up(void)
*****/
void	local_lookup(void)
{
  int i, j, cnt=0, num_of_nn_atoms=0;
  double dr=0.0;
  double dx, dy, dz;
  double rdx, rdy, rdz;
  double rxi, ryi, rzi;

  rdx = sys.Lx*0.5;
  rdy = sys.Ly*0.5;
  rdz = sys.Lz*0.5;

  for(i=0; i<sys.N; i++) {  /* "i" loop */
    rxi = sys.rx[i];
    ryi = sys.ry[i];
    rzi = sys.rz[i];

    num_of_nn_atoms = 0;

    for(j=0; j<sys.N; j++) {  /* "j" loop */

      if(i==j) continue;

      dx = rxi - sys.rx[j];
      dy = ryi - sys.ry[j];
      dz = rzi - sys.rz[j];

      /* cyclic boundary condition */

      if(dx > rdx) dx -= sys.Lx;
      if(dy > rdy) dy -= sys.Ly;
      if(dz > rdz) dz -= sys.Lz;
      if(dx < -rdx) dx += sys.Lx;
      if(dy < -rdy) dy += sys.Ly;
      if(dz < -rdz) dz += sys.Lz;

      dr = dx*dx + dy*dy + dz*dz; /* dr' = dr^2 */

      if(dr<local.Rc2) {
	local.lookup[cnt] = j; /* nn atom lookup list */
	cnt++;
	num_of_nn_atoms++;
      }

    }
    local.address[i+1] = cnt; /* store the address of the local.lookup */

    /* total number of nn atoms including the target atom */
    local.num_of_nn_atoms[i] = num_of_nn_atoms + 1; 
    if(local.num_of_nn_atoms[i] > local.maxN) {
      printf(">> the radius of the localized region = %f\n",local.Rc);
      printf(">> allowed num_of_nn_atoms (local.maxN) = %d\n",local.maxN);
      printf(">> local.num_of_nn_atoms[%d] = %d\n",i,local.num_of_nn_atoms[i]);
      printf(">>     adjust the local.maxN, local.Rc in [tb_init.c]\n");
      exit(0);
    }

  }

/*
  for(i=0;i<sys.N;i++) {
    printf("nn atoms of %d \n",i);
    printf("total number = %d\n",local.num_of_nn_atoms[i]);
    for(j=local.address[i];j<local.address[i+1];j++) {
      printf("%d ",local.lookup[j]);
    }
    printf("\n");
  }
  exit(0);
*/
}
