int calculate_glcm(int *image, char *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *glcm, int Ng);
int calculate_gldm(int *image, char *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *gldm, int Ng, int alpha);
int calculate_ngtdm(int *image, char *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *ngtdm, int Ng);
int calculate_glszm(int *image, signed char *mask, int Sx, int Sy, int Sz, int *angles, int Na, int *tempData, int Ng, int Ns);
int fill_glszm(int *tempData, double *glszm, int Ng, int maxRegion);