#include <Python.h>
#include "cmatrices.h"

int calculate_glcm(int *image, char *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *glcm, int Ng)
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

int calculate_gldm(int *image, char *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *gldm, int Ng, int alpha)
{
    /* Calculate GLDM: Count the number of voxels with gray level i, that have j dependent neighbours.
    *  A voxel is considered dependent if the absolute difference between the center voxel and the neighbour <= alpha
    */
    int gldm_idx_max = Ng * (Na * 2 + 1);
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
                    for (int d = 1; d >= -1; d-=2) // Iterate over both directions for each angle
	                {
                        for (int a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
                        {
                            // Check whether the neighbour index is not out of range (i.e. part of the image)
                            if (iz + d * angles[a * 3] >= 0 && iz + d * angles[a * 3] < Sz &&
                                iy + d * angles[a * 3 + 1] >= 0 && iy + d * angles[a * 3 + 1] < Sy &&
                                ix + d * angles[a * 3 + 2] >= 0 && ix + d * angles[a * 3 + 2] < Sx)
                            {
                                j = i + d * angles[a * 3 + 2] +
                                        d * angles[a * 3 + 1] * Sx +
                                        d * angles[a * 3] * Sy * Sx;
                                if (mask[j])  // Check whether neighbour voxel is part of the segmentation
                                {
                                    int diff = image[i] - image[j];
                                    if (diff < 0) diff *= -1;  // Get absolute difference
                                    if (diff <= alpha) dep++;
                                }
                            }
                        }
                    }
                    int gldm_idx = dep + (image[i]-1) * (Na * 2 + 1);
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

int calculate_ngtdm(int *image, char *mask, int Sx, int Sy, int Sz, int *angles, int Na, double *ngtdm, int Ng)
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
                    for (int d = 1; d >= -1; d-=2) // Iterate over both directions for each angle
	                {
                        for (int a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
                        {
                            // Check whether the neighbour index is not out of range (i.e. part of the image)
                            if (iz + d * angles[a * 3] >= 0 && iz + d * angles[a * 3] < Sz &&
                                iy + d * angles[a * 3 + 1] >= 0 && iy + d * angles[a * 3 + 1] < Sy &&
                                ix + d * angles[a * 3 + 2] >= 0 && ix + d * angles[a * 3 + 2] < Sx)
                            {
                                j = i + d * angles[a * 3 + 2] +
                                        d * angles[a * 3 + 1] * Sx +
                                        d * angles[a * 3] * Sy * Sx;
                                if (mask[j])  // Check whether neighbour voxel is part of the segmentation
                                {
                                   count++;
                                   sum += image[j];
                                }
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

int calculate_glszm(int *image, signed char *mask, int Sx, int Sy, int Sz, int *angles, int Na, int *tempData, int Ng, int Ns)
{
	int Ni = Sz * Sy * Sx;
	int maxRegion = 0;

	int regionCounter = 0;
	int max_region_idx = Ns * 2;
	int i = 0, j = 0, cur_idx = 0;
	for (int iz = 0; iz < Sz; iz++)
	{
		for (int iy = 0; iy < Sy; iy++)
		{
			for (int ix = 0; ix < Sx; ix++)
			{
				// Check if the current voxel is part of the segmentation and unprocessed
				if (mask[i] > 0)
				{
					// Store current gray level, needed for region growing and determines index in GLSZM.
					int gl = image[i];

					// Instantiate variable to hold region size at 1. Function region increases size for every found neighbour.
					// The initial size of 1 signifies current voxel
					int region = 0;

					// Start growing the region
					cur_idx = i;

					while (cur_idx > -1)
					{
						// Exclude current voxel to prevent reprocessing -> 0 is processed or not part of segmentation
						// -1 is marked for further processing
						mask[cur_idx] = 0;
						// Increment region size, as number of loops corresponds to number of voxels in current region
						region++;

						// Generate neighbours for current voxel
						int cur_z = (cur_idx / (Sx * Sy));
						int cur_y = (cur_idx % (Sx * Sy)) / Sx;
						int cur_x = (cur_idx % (Sx * Sy)) % Sx;
						for (int d = 1; d >= -1; d -= 2) // Iterate over both directions for each angle
						{
							for (int a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
							{
								// Check whether the neighbour index is not out of range (i.e. part of the image)
								if (cur_z + d * angles[a * 3] >= 0 && cur_z + d * angles[a * 3] < Sz &&
									cur_y + d * angles[a * 3 + 1] >= 0 && cur_y + d * angles[a * 3 + 1] < Sy &&
									cur_x + d * angles[a * 3 + 2] >= 0 && cur_x + d * angles[a * 3 + 2] < Sx)
								{
									j = cur_idx + d * angles[a * 3 + 2] +
										d * angles[a * 3 + 1] * Sx +
										d * angles[a * 3] * Sy * Sx;
									// Check whether neighbour voxel is part of the segmentation and unprocessed
									if (mask[j] && (image[j] == gl))
									{
										// Voxel belongs to current region, mark it for further processing
										mask[j] = -1;
									}
								}
							} // next a
						} // next d

						if (mask[j] == -1)
						{
						    cur_idx = j;
						}
						else
						{
						    int ref_idx = cur_idx;
						    cur_idx ++; // cur_idx doesn't have to be checked
						    // try to find the next voxel that has been marked for processing (-1)
                            // If none are found, current region is complete
                            // increment cur_idx until end of image or a voxel that has been marked is found
                            while (cur_idx < Ni && !(mask[cur_idx] == -1)) cur_idx++;
                            if (cur_idx == Ni) cur_idx = ref_idx; // no voxel found, check between ref_idx and i
                            {
                                while (cur_idx > i && !(mask[cur_idx] == -1)) cur_idx--;
                                if (cur_idx == i) cur_idx = -1;
                            }
                        }
					} // while region

					if (regionCounter >= max_region_idx) return -1; // index out of range

					if (region > maxRegion) maxRegion = region;

					tempData[(regionCounter * 2)] = gl;
					tempData[((regionCounter * 2) + 1)] = region;

					regionCounter ++;
				}
				// Increase index to point to next item. Only one index is needed as data is stored one dimensionally in memory.
				// This is in row-column-slice order.
				i++;
			}
		}
	}
	return maxRegion;
}

int fill_glszm(int *tempData, double *glszm, int Ng, int maxRegion)
{
    int i = 0;
	int glszm_idx_max = Ng * maxRegion;

    while(tempData[i * 2] > -1)
    {
        int glszm_idx = (tempData[i * 2] - 1) * maxRegion + tempData[i * 2 + 1] - 1;
        if (glszm_idx >= glszm_idx_max) return 0; // Index out of range

        glszm[glszm_idx]++;
        i++;
    }
    return 1;
}