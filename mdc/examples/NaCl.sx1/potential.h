
typedef	struct {
  /* ------------------------------------------------------------- */
  /* SX-1 potential parameter set: MEMO                            */
  /*                                                               */
  /*  function form:                                               */
  /*  u_(i,j) = Z_(i)*Z_(j)*e^(2)/r_(i,j) +                        */
  /*     f(b_(i)+b_(j))*Exp[(a_(i)+a_(j)-r_(i,j))/(b_(i)+b_(j))    */
  /*                                                               */
  /*  f0=1 [kcal/(A*mol)]                                          */ 
  /*  1 [kcal] = 4.1868*10^3  [J=N*m]                              */
  /*  1 [A] = 10^(-10) [m]                                         */
  /*                                                               */
  /*  for NaCl crystal system:                                     */
  /*  alpha[Na] = 1.260 [A]                                        */ 
  /*  alpha[Cl] = 1.950 [A]                                        */
  /*  beta[Na]  = 0.080 [A]                                        */
  /*  beta[Cl]  = 0.090 [A]                                        */
  /* ------------------------------------------------------------- */

  double *alpha;
  double *beta;
  double f0;
} potentitial_sx1_set;

/*------------------- declaration for the structures ----------------------*/
  potentitial_sx1_set sx1;

