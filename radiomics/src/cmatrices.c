#include <C:\Python27\include\Python.h>
#include "cmatrices.h"

int calculate_glcm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *glcm, int Ng)
{
    /* Calculate GLCM: Count the number of voxels with gray level i is neighboured by a voxel with gray level j in
    *  direction and distance specified by a. Returns an asymmetrical GLCM matrix for each angle/distance
    *  defined in angles.
    */
    int glcm_idx_max = Ng * Ng * Na;
    int i = 0, j = 0;
    for (int iz = 0; iz < Sz; iz ++) // Iterate over all voxels in a row - column - slice order
    {
        for (int iy = 0; iy < Sy; iy ++)
        {
            for (int ix = 0; ix < Sx; ix ++)
            {
                if (mask[i])  // Check if the current voxel is part of the segmentation
                {
                    for (int a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
                    {
                        // Check whether the neighbour index is not out of range (i.e. part of the image)
                        if (iz + angles[a * 3] >= 0 && iz + angles[a * 3] < Sz &&
                            iy + angles[a * 3 + 1] >= 0 && iy + angles[a * 3 + 1] < Sy &&
                            ix + angles[a * 3 + 2] >= 0 && ix + angles[a * 3 + 2] < Sx)
                        {
                            j = i + angles[a * 3 + 2] + angles[a * 3 + 1] * Sx + angles[a * 3] * Sy * Sx;
                            if (mask[j])  // Check whether neighbour voxel is part of the segmentation
                            {
                                int glcm_idx = a + (image[j]-1) * Na + (image[i]-1) * Na * Ng;
                                if (glcm_idx >= glcm_idx_max)
                                {
                                    return 0; // Index out of range
                                }
                                glcm[glcm_idx] ++;
                            }
                        }
                    }
                }
                // Increase index to point to next item. Only one index is needed as data is stored one dimensionally in memory.
                // This is in row-column-slice order.
                i++;
            }
        }
    }
    return 1;
}

int calculate_gldm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *gldm, int Ng, int Nd, int alpha)
{
    /* Calculate GLDM: Count the number of voxels with gray level i, that have j dependent neighbours.
    *  A voxel is considered dependent if the absolute difference between the center voxel and the neighbour <= alpha
    */
    int gldm_idx_max = Nd * Ng;
    int i = 0, j = 0;
    for (int iz = 0; iz < Sz; iz ++) // Iterate over all voxels in a row - column - slice order
    {
        for (int iy = 0; iy < Sy; iy ++)
        {
            for (int ix = 0; ix < Sx; ix ++)
            {
                if (mask[i])  // Check if the current voxel is part of the segmentation
                {
                    int dep = 0;
                    for (int a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
                    {
                        // Check whether the neighbour index is not out of range (i.e. part of the image)
                        if (iz + angles[a * 3] >= 0 && iz + angles[a * 3] < Sz &&
                            iy + angles[a * 3 + 1] >= 0 && iy + angles[a * 3 + 1] < Sy &&
                            ix + angles[a * 3 + 2] >= 0 && ix + angles[a * 3 + 2] < Sx)
                        {
                            j = i + angles[a * 3 + 2] + angles[a * 3 + 1] * Sx + angles[a * 3] * Sy * Sx;
                            if (mask[j])  // Check whether neighbour voxel is part of the segmentation
                            {
                                int diff = image[i] - image[j];
                                if (diff < 0) diff *= -1;  // Get absolute difference
                                if (diff <= alpha) dep++;
                            }
                        }
                    }
                    int gldm_idx = dep + (image[i]-1) * Nd;
                    if (gldm_idx >= gldm_idx_max) return 0; // Index out of range
                    gldm[gldm_idx] ++;
                }
                // Increase index to point to next item. Only one index is needed as data is stored one dimensionally in memory.
                // This is in row-column-slice order.
                i++;
            }
        }
    }
    return 1;
}

int calculate_ngtdm(int *image, int *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *ngtdm, int Ng)
{
    // Fill gray levels (empty slices gray levels are later deleted in python)
    for (int gl = 0; gl < Ng; gl++)
    {
        ngtdm[gl*3 + 2] = gl + 1;
    }

    /* Calculate matrix: for each gray level, element 0 describes the number of voxels with gray level i and
    *  element 1 describes the sum of all differences between voxels with gray level i and their neighbourhood
    */
    int ngtdm_idx_max = Ng * 3;
    int i = 0, j = 0;
    for (int iz = 0; iz < Sz; iz ++) // Iterate over all voxels in a row - column - slice order
    {
        for (int iy = 0; iy < Sy; iy ++)
        {
            for (int ix = 0; ix < Sx; ix ++)
            {
                if (mask[i])  // Check if the current voxel is part of the segmentation
                {
                    double count = 0;
                    double sum = 0;
                    for (int a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
                    {
                        // Check whether the neighbour index is not out of range (i.e. part of the image)
                        if (iz + angles[a * 3] >= 0 && iz + angles[a * 3] < Sz &&
                            iy + angles[a * 3 + 1] >= 0 && iy + angles[a * 3 + 1] < Sy &&
                            ix + angles[a * 3 + 2] >= 0 && ix + angles[a * 3 + 2] < Sx)
                        {
                            j = i + angles[a * 3 + 2] + angles[a * 3 + 1] * Sx + angles[a * 3] * Sy * Sx;
                            if (mask[j])  // Check whether neighbour voxel is part of the segmentation
                            {
                               count++;
                               sum += image[j];
                            }
                        }
                    }
                    double diff;
                    if (count == 0) { diff = 0; }
                    else { diff = (double)image[i] - sum  / count; }

                    if (diff < 0) diff *= -1;  // Get absolute difference

                    int ngtdm_idx = (image[i]-1) * 3;
                    if (ngtdm_idx >= ngtdm_idx_max) return 0; // Index out of range
                    ngtdm[ngtdm_idx]++;
                    ngtdm[ngtdm_idx + 1] += diff;
                }
                // Increase index to point to next item. Only one index is needed as data is stored one dimensionally in memory.
                // This is in row-column-slice order.
                i++;
            }
        }
    }
    return 1;
}