#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#include  "headers/md.h"

int  main(void)
{
  open_files();

  init_atomic_unit();
  set_atoms();
  set_param();

  dump_position_mesa(); /* dump initial position */

  init_bhs_pseudo_potential();
  /* check_bhs_dV_in_real_space(0, 0); */
  
  reciprocal_lattice_vector_generator();
  mkindex_fft3d();

  /* show_vectors(wv.nplw); */

  /*--- set Hamiltonian matrix of kinetic and pseudopotential ---*/
  clear_hmat_small();
  mk_Bessel(wv.nplw, wv.kgx, wv.kgy, wv.kgz);
  set_Vps_hmat(wv.nplwin, wv.kgx, wv.kgy, wv.kgz);
  set_Ek(wv.nplwin, wv.kgx, wv.kgy, wv.kgz);

  get_eigen();

  set_fermi();  /*--- set occupation number ---*/

  /*--- add hartree and correlation and exchange potential ---*/
  clear_hmat_small();
  set_Vps_hmat(wv.nplwin, wv.kgx, wv.kgy, wv.kgz);
  set_Ek(wv.nplwin, wv.kgx, wv.kgy, wv.kgz);

  get_chgdns();
  dump_chgdns();
  set_Vxc(wv.nplwin);
  set_Hartree(wv.nplwin);
  set_HVxc_hmat(wv.nplwin, wv.ngx, wv.ngy, wv.ngz);

  get_eigen();
  /*-------------------------------------------------------------*/
  
  init_plane_wave_coefficients();  

  init_ewald();
  ewald();

  Born_Oppenheimer();  /* approaching to the BO potential surface */
  set_force();

  /*================================================ MD start ===*/
  display(0);
  for(cpa.step=1; cpa.step<=cpa.max_step; cpa.step++){

    ewald();
    set_force();  /* F(t+dt), store F(t) for next_v()  */

    next_r();     /* Verlet's velocity form -> r(t+dt) */
    next_v();     /* Verlet's velocity form -> v(t+dt) */
    next_wf();    /* C(t+dt) */

    /* control temperature with velocity scaling */
    /* ctl_temp(ctl.step, ctl.temp); */

    get_energy();
    display(3);

    if(cp.residual_max > cpa.residual) {
      Born_Oppenheimer();   /* approaching to the BO potential surface */
      set_force();
    }

    /*--- dump data ---*/
    dump_eigen();
    dump_position();
    dump_energy();
    dump_force();
    dump_velocity();
    dump_temp();
    dump_position_mesa();
    /*-----------------*/
    
  } /*======================================== End of MD Loop ===*/  

  close_files();  
  return 0;
}
