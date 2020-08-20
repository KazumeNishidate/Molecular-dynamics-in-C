#include  <stdlib.h>
#include  <stdio.h>
#include  "md.h"
#include  "prototypes.h"

void calc_R_vs_E(void) {
  int num_divisions;
  double lowest_distance, max_distance;
  double delta_distance;

  /*---------------------------------*/
  lowest_distance = 0.8;  /* [A] */
  max_distance = 2.8;
  num_divisions = 100;
  /*---------------------------------*/

  get_control_param();  /* calculational control parameters: control.c */

  sys.Ax = lowest_distance;
  sys.Ay = lowest_distance;
  sys.Az = lowest_distance;

  delta_distance = (max_distance-lowest_distance)/(double)num_divisions;

  init_mem();           /* dynamic memory allocation        */
  set_potential();      /* potential setting: control.c     */
  identify_ion();       /* identify the ions: control.c     */

  while(sys.Ax < max_distance) {
    sys.Ax += delta_distance;
    sys.Ay += delta_distance;
    sys.Az += delta_distance;

    /*========= TBMD ==============*/
    set_roc();            /* initial locations: control.c     */
    clear_foc();          /* clear the forces                 */
    real_space();         /*  real space sum for force and potential     */
    jacobi_trans();             /* diagonalization */
    calc_elec_E();
    /*=============================*/

    printf("%f   %f\n",sys.Ax,sys.pot/(double)sys.N);

  }
  exit(0);   

}

void show_Htb(void) {
  int i, j;

  printf("\n");
  printf("--- Tight-Binding Hamiltonian --- \n");
  for(i=0;i<sys.N*4;i++)
    {
      for(j=0;j<sys.N*4;j++)
        {
          printf("%8.4f",htb.mat1[sys.N*4*i+j]);
        }
      printf("\n");
    }
  getchar();

  printf("--- force matrix for x direction --- \n");
  for(i=0;i<sys.N*4;i++)
    {
      for(j=0;j<sys.N*4;j++)
        {
          printf("%8.4f",htb.mat2x[sys.N*4*i+j]);
        }
      printf("\n");
    }
  getchar();
}

