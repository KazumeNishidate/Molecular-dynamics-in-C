/*********************************************************************/
/*  Silicon (diamond structure) with Stillinger Weber potential      */
/*********************************************************************/

typedef struct{
  double ep;   /* energy unit    */
  double sig;  /* length unit    */
  double A;    /* parameter set  */
  double B;
  double a;
  double lam;
  double gam;
  double beta;

  double f1x, f1y, f1z;
  double f2x, f2y, f2z;
  double f3x, f3y, f3z;
  double f4x, f4y, f4z;
  double f5x, f5y, f5z;
  double f6x, f6y, f6z;
  double f7x, f7y, f7z;
  double f8x, f8y, f8z;
  double f9x, f9y, f9z;
  double f1_a, f2_a, f5_a;
  double f1_b, f2_b, f5_b;

} potential_sw_set;

/*------------------- declaration for the structures ----------------------*/
  potential_sw_set sw;
