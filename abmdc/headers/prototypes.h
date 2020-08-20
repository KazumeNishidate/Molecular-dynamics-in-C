/*----------- prototype declarations -------------*/

/* main.c */
int   main(void);

/* control.c */
void   init_bhs_pseudo_potential(void);
void   set_param(void);
void   set_atoms(void);
void   set_rocs(void);

/* bhs_init.c */
void   check_bhs_dV_in_real_space(int atom, int eru);
void   init_bhs_c2a(int atom);

/* plane_wave.c */
void   orthogonality_check(void);
void   mkindex_fft3d(void);
void   reciprocal_lattice_vector_generator(void);
void   init_plane_wave_coefficients(void);
void   randomize(void);
void   normalize(double *, double *);
void   show_vectors(int npl);
void   check_array_size(int cnt, int cnti);

/* ext.c */
double   erfcc(double x);

/* msd.c */
void calc_msd(void);

/* init.c */
void   init_atomic_unit(void);
void   init_mem(void);
void   init_velocity(void);

/* pseudo.c */
void   set_Vps(int npl, double *kx, double *ky, double *kz);
void   get_HF_force(int n_pl, double *kx, double *ky, double *kz);
void   clear_hmat_small(void);
void   set_Ek(int npw, double *kx, double *ky, double *kz);
void   set_Vps_hmat(int npl, double *kx, double *ky, double *kz);
void   mk_Bessel(int n_pl, double *kx, double *ky, double *kz);
void   calc_core_Vps(int type, int nq1, int nq2, 
		     double *kx, double *ky, double *kz);
void   calc_delta_Vps(int type, int nq1, int nq2,
		      double *kx, double *ky, double *kz);
void   calc_delta_Vps_zero(int type, int nq1, int nq2,
			   double *kx, double *ky, double *kz);

/* eigen.c */
void  get_eigen(void);
void  check_ckg_small(void);
void  househ(void);
void  bisect(void);
void  invitr(void);
void  invitr_vh(void);
void  check_hamiltonian_small(void);
void  check_eigen_values_small(void);

/* fermi.c */
void  set_fermi(void);
void  get_chgdns(void);

/* fft.c */
void  fft3d(int isign);
void  check_fft(void);

/* ewald.c */
void  init_ewald(void);
double find_min3(double x, double y, double z);
void  ewald(void);
void  clear_fcc(void);
void  ewald_I(void);
void  ewald_II(void);
void  ewald_III(void);

/* xc.c */
void  set_HVxc_hmat(int npw, int *nx, int *ny, int *nz);
void  set_Vxc(int npw);
void  calc_Vxc(void);
void  set_Hartree(int npw);
void  set_VhxcC(int npw, int *nx, int *ny, int *nz);

/* energy.c */
void  get_energy(void);
void  calc_kinetic(void);

/* quench.c */
void  Born_Oppenheimer(void);
void  set_Legendre();
void  quench_wv(void);
void  mk_HxC(void);
void  Gram_Schmidt(void);
void  sort_eigen();
void  show_eigen();
double   inner_productR(int ib, int jb);
double   inner_productI(int ib, int jb);

/* next.c */
void  next_wf(void);
void  next_r(void);  /* Verlet */
void  next_v(void);  /* Verlet */
void  Ryckaert(void);
void  get_next_wf(void);
void  set_force(void);

/* pt.c */
void  ctl_temp(int step, double Tset);
void  calc_press(void);

/* monitor.c */
void  open_files(void);
void  close_files(void);
void  dump_eigen(void);
void  dump_chgdns(void);
void  dump_position(void);
void  dump_energy(void);
void  dump_force(void);
void  dump_velocity(void);
void  dump_temp(void);
void  dump_pres(void);
void  dump_position_mesa(void);
void  dump_msd(void);
void  display_md_start(void);
void  display_bo_start(void);
void  display_bo_end(void);
void  display(void);
