int calculate_glcm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *glcm, int Ng);
int calculate_gldm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *gldm, int Ng, int Nd, int alpha);
int calculate_ngtdm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *ngtdm, int Ng);