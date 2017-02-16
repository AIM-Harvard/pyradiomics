int calculate_glcm(int *image, char *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *glcm, int Ng);
int calculate_glszm(int *image, signed char *mask, int Sx, int Sy, int Sz, int *angles, int Na, int *tempData, int Ng, int Ns);
int fill_glszm(int *tempData, double *glszm, int Ng, int maxRegion);
int calculate_glrlm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *glrlm, int Ng, int Nr);
int run_diagonal(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *glrlm, int glrlm_idx_max, int Nr, int *jd, int a);
