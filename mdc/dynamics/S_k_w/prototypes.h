/* main.c */
int	main(void);
void	calc_Skw(void);

/* files.c */
void	openfiles(void);
void 	close_files(void);
void   record_data(void);

/* fft1.c */
void fft(double* x, double* y, int step, double flag);
void hanning_window(double *x, int n);

/* control.c */
void  init_param(void);
void  init_mem(void);
