/*----------- prototype declarations -------------*/

/* control.c */
void   read_poscar(void);
void   get_control_param(void);
void   identify_ion(void);
void   set_roc(void);
void   set_potential(void);

/* ext.c */
double   nrand(double temp, double bunsan);

/* files.c */
void   open_files(void);
void   close_files(void);
void   print_to_file(void);
void   display4(void);
void   r_position(void);
void   r_velocities(void);

/* init.c */
void   init_mem(void);
void   set_vel(void);
void   clear_foc(void);
void   init(void);

/* pt.c */
void   calc_press(void);
void   control_press(int d_step);
void   control_temp(int d_step, double temp);

/* real.c */
void   real_space(void);
void   sort_eigenvalue(void);
void   calc_elec_E(void);

/* jacobi.c */
void   jacobi_trans(void);

/* rv.c */
void   calc_foc(void);
void   next_rv_verlet(void);
void   next_r(void);
void   next_v(void);
void   calc_kin(void);
void   moment_correction(void);
void   next_rv_gear(void);

/* msd.c */
void   calc_msd(void);

/* debug.c */
void show_Htb(void);
void calc_R_vs_E(void);
