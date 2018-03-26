int calculate_glcm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *glcm, int Ng);
int calculate_glszm(int *image, char *mask, int *size, int *strides, int *angles, int Na, int *tempData, int Ng, int Ns);
int fill_glszm(int *tempData, double *glszm, int Ng, int maxRegion);
int calculate_glrlm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *glrlm, int Ng, int Nr);
int run_diagonal(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *glrlm, int glrlm_idx_max, int Nr, int *jd, int a);
int calculate_ngtdm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *ngtdm, int Ng);
int calculate_gldm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *gldm, int Ng, int alpha);
int generate_angles(int *size, int *distances, int n_dim, int n_dist, char bidirectional, int force2Ddim, int **angles, int *n_a);
