int calculate_glcm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *glcm, int Ng);
int calculate_gldm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *gldm, int Ng, int Nd, int alpha);
int calculate_ngtdm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *ngtdm, int Ng);
int calculate_glszm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *glszm, int Ng, int Ns);
void grow_region(int *image, int *mask, int Sx, int Sy, int Sz, int x, int y, int z, int *angles, int Na, int gl, int *regionptr);