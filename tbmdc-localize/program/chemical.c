#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  "md.h"
#include  "prototypes.h"

/*===================================================================*/
#define PRECISION 10.0   /* precision of the N_e calculation [SWNT:100.0] */
#define RATE_OF_CHANGE_MYU 3.0   /* rate of change for ch.myu [SWNT:10.0] */
#define DISP_INFO

void set_chemical_potential( void ) {

  int cnt = 0, N_e, iN_e_calc = 0;
  double delta_myu = 1.0*htb.eV2E;  /* = 1.0 [eV] */
  /* double old_myu = ch.myu; */
  double pre_myu = 0.0, pre_myu2 = 0.0;
  double dN_e_calc = 0.0, diff;

  /* N_e = number of electrons [particle_number * 2] * PRECISION   */
  N_e = (int)fermi.electrons*PRECISION;   /* spin degeneracy is considered */
  dN_e_calc = calc_number_of_electron();
  iN_e_calc = (int)(dN_e_calc);

  diff = fabs( (double)(dN_e_calc - N_e)/(double)PRECISION );

  while(diff<1.0){
    delta_myu /= RATE_OF_CHANGE_MYU;
    diff *= 10.0;
  }
  delta_myu /= RATE_OF_CHANGE_MYU;

  while((N_e!=iN_e_calc)&&cnt<1000){

#ifdef DISP_INFO
    printf(">> %d\n",cnt);
    printf("   (ch.myu : delta_myu  ) = (%f : %f) [eV]\n",
	    ch.myu/htb.eV2E,delta_myu/htb.eV2E );
    printf("   (N_e    : [N_e_calc] ) = (%d : [%d : %f])\n", N_e, iN_e_calc,dN_e_calc );
#endif

    pre_myu2 = pre_myu;
    pre_myu  = ch.myu;

    if(N_e < iN_e_calc) ch.myu -= delta_myu;
    if(N_e > iN_e_calc) ch.myu += delta_myu;

    if(fabs(ch.myu-pre_myu2) < 10e-17){
      /* if(ch.myu == pre_myu2){ */
      delta_myu /= RATE_OF_CHANGE_MYU;
      ch.myu  = pre_myu;
      if( ch.myu > pre_myu2 ) ch.myu -= delta_myu;
      else                    ch.myu += delta_myu;
    }
    dN_e_calc = calc_number_of_electron();
    iN_e_calc = (int)dN_e_calc;
    cnt++;
  }

  if( cnt >= 1000 ){
    printf(">>> adjusting procedure is aborted\n");
    /* printf(">>> reset the chemical potential %f to %f\n",ch.myu,old_myu); */
    /* ch.myu = old_myu; */
  }

}

double calc_number_of_electron( void ) {
  int target_atom;
  double dN_e_calc = 0.0;

  get_chebyshev_coeff();
  for(target_atom = 0; target_atom < sys.N; target_atom++) {
    hamiltonian(target_atom);
    calc_Fermi_E(target_atom);
    dN_e_calc += (fermi.number_of_electrons*PRECISION);
  }

  return dN_e_calc;
}

/*===================================================================*/

