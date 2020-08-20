
typedef struct {
  /* --------------------------------------------------------------- */
  /*   Soft-core-potential                                           */
  /* --------------------------------------------------------------- */
  int nn, soft;
  double epsilon;
  double sigma_ca;
  double sigma_f;

} soft_core_potential_set;

soft_core_potential_set soft;
