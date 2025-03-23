int calculate_glcm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, double *glcm, int Ng);
int calculate_glszm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, int *tempData, int Ng, int Ns, int Nvox);
int fill_glszm(int *tempData, double *glszm, int Ng, int maxRegion);
int calculate_glrlm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, double *glrlm, int Ng, int Nr);
int calculate_ngtdm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, double *ngtdm, int Ng);
int calculate_gldm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, double *gldm, int Ng, int alpha);
int get_angle_count(int *size, int *distances, int Nd, int n_dist, char bidirectional, int force2Ddim);
int build_angles(int *size, int *distances, int Nd, int n_dist, int force2Ddim, int Na, int *angles);
