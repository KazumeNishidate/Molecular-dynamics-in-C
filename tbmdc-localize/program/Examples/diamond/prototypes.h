/*----------- prototype declarations -------------*/

/* main.c */
int   main(void);
void   newton(void);

/* control.c */
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

/* tb.c */
void   fill_htb1(int local_N, int i, int j, double dr, double s_R, 
		 double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_sssig_r, double v_spsig_r,  
		 double v_ppsig_r, double v_pppi_r);
void   fill_htb2(int local_N, int i, int j, double dr, double s_R, 
		 double s_dash_R, double eru, double emu, double enu,
		 double eru2, double emu2, double enu2, 
		 double v_sssig_r, double v_spsig_r,  
		 double v_ppsig_r, double v_pppi_r);

void   diagonal(int local_N, int i);

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

/* chebyshev.c */
void   get_chebyshev_coeff(void);

/* fermi.c */
void   calc_Fermi_F(int target_atom);
void   calc_Fermi_E(int target_atom);
void   check_H_chebyshev(void);

/* chemical.c */
void   set_chemical_potential( void );
double calc_number_of_electron( void );

/* tb_init.c */
void   init_mem_chebyshev_coeff(void);

/* lookup.c */
void local_lookup(void);

/* hamiltonian.c */
void	hamiltonian(int atom);

/* repulsive.c */
void	repulsive(void);

/* ./xsrc/network.c */
void   network(void);

/* ./xsrc/xdmain.c */
void   xd(void);
void   open_xd(void);

/* ./xsrc/xdicon.c */
void   xd_icon_set(void);
void   xd_icon_clear(void);
void   xd_icon_heat(void);
void   xd_icon_cool(void);
void   xd_icon_xhist(void);
void   xd_icon_xunit(void);
void   xd_icon_xnet(void);
void   xd_icon_xmsd(void);

/* ./xsrc/xdinit.c */
void   xd_init_set(void);

/* ./xsrc/xhist.c */
void   open_xhist(void);

/* ./xsrc/xmsd.c */
void   open_xmsd(void);

/* ./xsrc/xnet.c */
void   open_xnet(void);
void   xnet_drag(void);

/* ./xsrc/xunit.c */
void   open_xunit(void);
void   xunit_drag(void);

/* xmycolor.c */
unsigned long MyColor();

