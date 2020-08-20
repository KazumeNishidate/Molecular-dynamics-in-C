
typedef struct {
  /* --------------------------------------------------------------- */
  /*   Born-Mayer-Huggins potential with Tosi-Fumi parameter set     */
  /*   see set_potential_hm() in "control.c" and references          */
  /* --------------------------------------------------------------- */
  double b;
  double c_na_na;
  double c_na_cl;
  double c_cl_cl;

  double sigma_na;
  double sigma_cl;

  double rho;

  double C_na_na;
  double C_na_cl;
  double C_cl_cl;

  double D_na_na;
  double D_na_cl;
  double D_cl_cl;

} huggins_mayer_potential_set;

/*------------------- declaration for the structures ----------------------*/
  huggins_mayer_potential_set hm;


