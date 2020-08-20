/*----------- prototype declarations -------------*/

/* main.c */
int   main(void);
void   newton(void);

/* control.c */
void   read_poscar(void);
void   get_control_param(void);
void   unit_converter(void);
void   identify_ion(void);
void   set_loc(void);
void   set_potential(void);
void   mk_table(void);

/* ext.c */
double   nrand(double temp, double bunsan);
double   erfcc(double x);

/* files.c */
void   open_files(void);
void   close_files(void);
void   print_to_file(void);
void   display1(void);
void   display2(void);
void   no_display(void);
void   r_position(void);
void   r_velocities(void);
void   r_macmol(void);

/* init.c */
void   init_mem(void);
void   set_vel(void);
void   clear_foc(void);
void   init(void);
void   calc_alpha(void);          /* alpha estimation */
void   reciprocal_space3a(void);  /* alpha estimation */

/* pt.c */
void   calc_press(void);
void   control_press(int d_step);
void   control_temp(int d_step, double temp);

/* rcprcl.c */
void   reciprocal_space(void);
void   reciprocal_space3(void);

/* real.c */
void   real_space(void);
void   real_vdW(void);

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

