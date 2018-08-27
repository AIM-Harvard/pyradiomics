#include <stdlib.h>
#include "cmatrices.h"

int calculate_glcm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *glcm, int Ng)
{
  /* Calculate GLCM: Count the number of voxels with gray level i is neighboured by a voxel with gray level j in
  *  direction and distance specified by a. Returns an asymmetrical GLCM matrix for each angle/distance
  *  defined in angles.
  */
  int glcm_idx_max = Ng * Ng * Na;
  int i = 0, j = 0;
  int iz, iy, ix;
  int a, glcm_idx;
  for (iz = 0; iz < size[0]; iz ++) // Iterate over all voxels in a row - column - slice order
  {
    for (iy = 0; iy < size[1]; iy ++)
    {
      for (ix = 0; ix < size[2]; ix ++)
      {
        if (mask[i])  // Check if the current voxel is part of the segmentation
        {
          for (a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
          {
            // Check whether the neighbour index is not out of range (i.e. part of the image)
            if (iz + angles[a * 3] >= 0 && iz + angles[a * 3] < size[0] &&
                iy + angles[a * 3 + 1] >= 0 && iy + angles[a * 3 + 1] < size[1] &&
                ix + angles[a * 3 + 2] >= 0 && ix + angles[a * 3 + 2] < size[2])
            {
              j = i + angles[a * 3] * strides[0] +
                  angles[a * 3 + 1] * strides[1] +
                  angles[a * 3 + 2] * strides[2];
              if (mask[j])  // Check whether neighbour voxel is part of the segmentation
              {
                glcm_idx = a + (image[j]-1) * Na + (image[i]-1) * Na * Ng;
                if (glcm_idx >= glcm_idx_max) return 0; // Index out of range
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

int calculate_glszm(int *image, char *mask, int *size, int *strides, int *angles, int Na, int *tempData, int Ng, int Ns)
{
  int Ni = size[0] * size[1] * size[2];
  int maxSize = 0;

  int *regionStack;
  int stackTop = -1;

  int regionCounter = 0;
  int max_region_idx = Ns * 2;

  int gl, region;
  int i, a;
  int cur_idx, cur_x, cur_y, cur_z, j;

  regionStack = (int *)calloc(Ns, sizeof(int));

  for (i = 0; i < Ni; i++)  // Loop over all the voxels in the image
  {
    // Check if the current voxel is part of the segmentation and unprocessed
    if (mask[i])
    {
      // Store current gray level, needed for region growing and determines index in GLSZM.
      gl = image[i];

      // Instantiate variable to hold region size at 0. Region increases for every found voxel.
      region = 0;

      // Start growing the region
      regionStack[++stackTop] = i; // Add the current voxel to the stack as 'starting point'
      mask[i] = 0;  // Mark current voxel as 'processed'

      while (stackTop > -1)
      {
        cur_idx = regionStack[stackTop--];  // Get the next voxel to process, on first iteration, this equals i

        // Increment region size, as number of loops corresponds to number of voxels in current region
        region++;

        // Generate neighbours for current voxel
        cur_z = (cur_idx / strides[0]);
        cur_y = (cur_idx % strides[0]) / strides[1];
        cur_x = (cur_idx % strides[0]) % strides[1];
        for (a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
        {
          // Check whether the neighbour index is not out of range (i.e. part of the image)
          if (cur_z + angles[a * 3] >= 0 && cur_z + angles[a * 3] < size[0] &&
            cur_y + angles[a * 3 + 1] >= 0 && cur_y + angles[a * 3 + 1] < size[1] &&
            cur_x + angles[a * 3 + 2] >= 0 && cur_x + angles[a * 3 + 2] < size[2])
          {
            j = cur_idx + angles[a * 3] * strides[0] +
                          angles[a * 3 + 1] * strides[1] +
                          angles[a * 3 + 2] * strides[2];
            // Check whether neighbour voxel is part of the current region and unprocessed
            if (mask[j] && (image[j] == gl))
            {
              // Push the voxel index to the stack for further processing
              regionStack[++stackTop] = j;
              // Voxel belongs to current region, mark it as 'processed'
              mask[j] = 0;
            }
          }
        } // next a
      }  // while (stackTop > -1)

      if (regionCounter >= max_region_idx)
      {
        free(regionStack);
        return -1; // index out of range
      }
      if (region > maxSize) maxSize = region;

      tempData[(regionCounter * 2)] = gl;
      tempData[((regionCounter * 2) + 1)] = region;

      regionCounter++;
    }
  }
  free(regionStack);

  if (regionCounter >= max_region_idx) return -1; // index out of range
  tempData[(regionCounter * 2)] = -1; // Set the first element after last region to -1 to stop the loop in fill_glszm
  return maxSize;
}

int fill_glszm(int *tempData, double *glszm, int Ng, int maxRegion)
{
  int i = 0;
  int glszm_idx_max = Ng * maxRegion;
  int glszm_idx;

  while(tempData[i * 2] > -1)
  {
    glszm_idx = (tempData[i * 2] - 1) * maxRegion + tempData[i * 2 + 1] - 1;
    if (glszm_idx >= glszm_idx_max) return 0; // Index out of range

    glszm[glszm_idx]++;
    i++;
  }
  return 1;
}

int calculate_glrlm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *glrlm, int Ng, int Nr)
{
  int a, da, nd;
  int mDims[3], sDims[2];
  int d1, d2;
  int jd[3];
  int glrlm_idx_max = Ng * Nr * Na;
  int runVal;
  char multiElement;

  for (a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
  {
    multiElement = 0;
    nd = 0;
    for (da = 0; da < 3; da++)
    {
      if (angles[a * 3 + da] == 0)
      {
        sDims[da - nd] = da;
      }
      else
      {
        mDims[nd] = da;
        nd++;
      }
    }
    if (nd == 1)  // Moving in 1 dimension
    {
      // Iterate over all start voxels: All voxels, where moving dim == 0 or max (i.e. all combinations of voxels
      // in static dimensions)
      for (d1 = 0; d1 < size[sDims[0]]; d1++)
      {
        for (d2 = 0; d2 < size[sDims[1]]; d2++)
        {
          jd[sDims[0]] = d1;
          jd[sDims[1]] = d2;
          if (angles[a * 3 + mDims[0]]< 0)
          {
            jd[mDims[0]] = size[mDims[0]] - 1; // Moving dim is negative, set start at maximum index
          }
          else
          {
            jd[mDims[0]] = 0; // Moving dim is positive, set start at 0
          }
          // run_diagonal returns -1 if error, otherwise number of elements encountered
          runVal = run_diagonal(image, mask, size, strides, angles, Na, glrlm, glrlm_idx_max, Nr, jd, a);
          if (runVal < 0) return runVal; // index out of range, 'raise' error.
          if (!multiElement && runVal > 1) multiElement = 1; // multiple elements found
        }
      }
    }
    else if (nd == 2) // Moving in 2 dimensions
    {
      for (d1 = 0; d1 < size[sDims[0]]; d1++) // iterate over all voxels in the static dimension
      {
        for (d2 = 0; d2 < size[mDims[0]]; d2++) // Iterate over moving dimension 1, treating it as a 'static'
        {
          jd[sDims[0]] = d1;
          jd[mDims[0]] = d2;
          if (angles[a * 3 + mDims[1]] < 0)
          {
            jd[mDims[1]] = size[mDims[1]] - 1; // Moving dim is negative, set start at maximum index
          }
          else
          {
            jd[mDims[1]] = 0; // Moving dim is positive, set start at 0
          }
          // run_diagonal returns -1 if error, otherwise number of elements encountered
          runVal = run_diagonal(image, mask, size, strides, angles, Na, glrlm, glrlm_idx_max, Nr, jd, a);
          if (runVal < 0) return runVal; // index out of range, 'raise' error.
          if (!multiElement && runVal > 1) multiElement = 1; // multiple elements found
        }
        for (d2 = 1; d2 < size[mDims[1]]; d2++) // Iterate over other moving dimension, treating it as a 'static'
        {
          jd[sDims[0]] = d1;
          jd[mDims[1]] = d2;
          // Prevent calculating the true diagonal twice, which is either at index 0, or max, depending on the
          // angle of the the second moving dimension: iterate 1 to max if angle is positive, 0 to max -1 if
          // negative. The ignored index is handled in the previous loop.
          if (angles[a * 3 + mDims[1]] < 0) jd[mDims[1]] --;
          if (angles[a * 3 + mDims[0]] < 0)
          {
            jd[mDims[0]] = size[mDims[0]] - 1; // Moving dim is negative, set start at maximum index
          }
          else
          {
            jd[mDims[0]] = 0; // Moving dim is positive, set start at 0
          }
          // run_diagonal returns -1 if error, otherwise number of elements encountered
          runVal = run_diagonal(image, mask, size, strides, angles, Na, glrlm, glrlm_idx_max, Nr, jd, a);
          if (runVal < 0) return runVal; // index out of range, 'raise' error.
          if (!multiElement && runVal > 1) multiElement = 1; // multiple elements found
        }
      }
    }
    else if (nd == 3)
    {
      for (d1 = 0; d1 < size[mDims[0]]; d1++)
      {
        for (d2 = 0; d2 < size[mDims[1]]; d2++)
        {
          jd[mDims[0]] = d1;
          jd[mDims[1]] = d2;
          if (angles[a * 3 + mDims[2]] < 0)
          {
            jd[mDims[2]] = size[mDims[2]] - 1; // Moving dim is negative, set start at maximum index
          }
          else
          {
            jd[mDims[2]] = 0; // Moving dim is positive, set start at 0
          }
          // run_diagonal returns -1 if error, otherwise number of elements encountered
          runVal = run_diagonal(image, mask, size, strides, angles, Na, glrlm, glrlm_idx_max, Nr, jd, a);
          if (runVal < 0) return runVal; // index out of range, 'raise' error.
          if (!multiElement && runVal > 1) multiElement = 1; // multiple elements found
        }
      }
      for (d1 = 0; d1 < size[mDims[1]]; d1++)
      {
        for (d2 = 1; d2 < size[mDims[2]]; d2++)
        {
          jd[mDims[1]] = d1;
          jd[mDims[2]] = d2;
          // Prevent calculating the true diagonal twice, which is either at index 0, or max, depending on the
          // angle of the the second moving dimension: iterate 1 to max if angle is positive, 0 to max -1 if
          // negative. The ignored index is handled in the previous loop.
          if (angles[a * 3 + mDims[2]] < 0) jd[mDims[2]] --;
          if (angles[a * 3 + mDims[0]] < 0)
          {
            jd[mDims[0]] = size[mDims[0]] - 1; // Moving dim is negative, set start at maximum index
          }
          else
          {
            jd[mDims[0]] = 0; // Moving dim is positive, set start at 0
          }
          // run_diagonal returns -1 if error, otherwise number of elements encountered
          runVal = run_diagonal(image, mask, size, strides, angles, Na, glrlm, glrlm_idx_max, Nr, jd, a);
          if (runVal < 0) return runVal; // index out of range, 'raise' error.
          if (!multiElement && runVal > 1) multiElement = 1; // multiple elements found
        }
      }
      for (d1 = 1; d1 < size[mDims[2]]; d1++)
      {
        for (d2 = 1; d2 < size[mDims[0]]; d2++)
        {
          jd[mDims[2]] = d1;
          jd[mDims[0]] = d2;
          // Prevent calculating the true diagonal twice, which is either at index 0, or max, depending on the
          // angle of the the second moving dimension: iterate 1 to max if angle is positive, 0 to max -1 if
          // negative. The ignored index is handled in the previous loops. This is done for both 'statics'.
          if (angles[a * 3 + mDims[2]] < 0) jd[mDims[2]] --;
          if (angles[a * 3 + mDims[0]] < 0) jd[mDims[0]] --;
          if (angles[a * 3 + mDims[1]] < 0)
          {
            jd[mDims[1]] = size[mDims[1]] - 1; // Moving dim is negative, set start at maximum index
          }
          else
          {
            jd[mDims[1]] = 0; // Moving dim is positive, set start at 0
          }
          // run_diagonal returns -1 if error, otherwise number of elements encountered
          runVal = run_diagonal(image, mask, size, strides, angles, Na, glrlm, glrlm_idx_max, Nr, jd, a);
          if (runVal < 0) return runVal; // index out of range, 'raise' error.
          if (!multiElement && runVal > 1) multiElement = 1; // multiple elements found
        }
      }
    } // switch No of dimensions
    if (!multiElement) // Segmentation is 2D for this angle, remove it.
    {
      // set elements at runlength index 0 to 0 for all gray levels in this angle, other runlength indices are
      // already 0.
      for (d1 = 0; d1 < Ng; d1++)
      {
        glrlm[a + d1 * Na * Nr] = 0;
      }
    }
  } // next a
  return 1;
}

int run_diagonal(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *glrlm, int glrlm_idx_max, int Nr, int *jd, int a)
{
  int j, gl, rl;
  int glrlm_idx;
  int elements = 0;

  gl = -1;
  rl = 0;

  j = jd[0] * strides[0] + jd[1] * strides[1] + jd[2] * strides[2];
  while (jd[0] >= 0 && jd[0] < size[0] &&
        jd[1] >= 0 && jd[1] < size[1] &&
        jd[2] >= 0 && jd[2] < size[2])
  {
    j = jd[0] * strides[0] + jd[1] * strides[1] + jd[2] * strides[2];
    if (mask[j])
    {
        elements++; // Count the number of segmented voxels in this run
      if (gl == image[j])
      {
        rl++;
      }
      else if (gl == -1)
      {
        gl = image[j];
        rl = 0;
      }
      else
      {
        glrlm_idx = a + rl * Na + (gl - 1) * Na * Nr;
        if (glrlm_idx >= glrlm_idx_max) return -1;
        glrlm[glrlm_idx]++;

        gl = image[j];
        rl = 0;
      }
    }
    else if (gl > -1)
    {
      glrlm_idx = a + rl * Na + (gl - 1) * Na * Nr;
      if (glrlm_idx >= glrlm_idx_max) return -1;
      glrlm[glrlm_idx]++;
      gl = -1;
    }
    jd[0] += angles[a * 3];
    jd[1] += angles[a * 3 + 1];
    jd[2] += angles[a * 3 + 2];
  }
  if (gl > -1)
  {
    glrlm_idx = a + rl * Na + (gl - 1) * Na * Nr;
    if (glrlm_idx >= glrlm_idx_max) return -1;
    glrlm[glrlm_idx]++;
  }
  return elements;
}

int calculate_ngtdm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *ngtdm, int Ng)
{
  int gl;
  int ngtdm_idx_max = Ng * 3;
  int i = 0, j = 0;
  int iz, iy, ix;
  int offset[3];
  double count, sum, diff;
  int a;
  int ngtdm_idx;

  // Fill gray levels (empty slices gray levels are later deleted in python)
  for (gl = 0; gl < Ng; gl++)
  {
      ngtdm[gl*3 + 2] = gl + 1;
  }
  /* Calculate matrix: for each gray level, element 0 describes the number of voxels with gray level i and
  *  element 1 describes the sum of all differences between voxels with gray level i and their neighbourhood
  */
  for (iz = 0; iz < size[0]; iz ++) // Iterate over all voxels in a row - column - slice order
  {
    for (iy = 0; iy < size[1]; iy ++)
    {
      for (ix = 0; ix < size[2]; ix ++)
      {
        if (mask[i])  // Check if the current voxel is part of the segmentation
        {
          count = 0;
          sum = 0;
          for (a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
          {
            offset[0] = angles[a * 3];
            offset[1] = angles[a * 3 + 1];
            offset[2] = angles[a * 3 + 2];

            // Check whether the neighbour index is not out of range (i.e. part of the image)
            if (iz + offset[0] >= 0 && iz + offset[0] < size[0] &&
                iy + offset[1] >= 0 && iy + offset[1] < size[1] &&
                ix + offset[2] >= 0 && ix + offset[2] < size[2])
            {
              j = i + offset[0] * strides[0] +
                      offset[1] * strides[1] +
                      offset[2] * strides[2];
              if (mask[j])  // Check whether neighbour voxel is part of the segmentation
              {
               count++;
               sum += image[j];
              }
            }
          }
          if (count == 0) { diff = 0; }
          else { diff = (double)image[i] - sum  / count; }

          if (diff < 0) diff *= -1;  // Get absolute difference

          ngtdm_idx = (image[i]-1) * 3;
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

int calculate_gldm(int *image, char *mask, int *size, int *strides, int *angles, int Na, double *gldm, int Ng, int alpha)
{
  /* Calculate GLDM: Count the number of voxels with gray level i, that have j dependent neighbours.
  *  A voxel is considered dependent if the absolute difference between the center voxel and the neighbour <= alpha
  */
  int gldm_idx_max = Ng * (Na * 2 + 1);
  int i = 0, j = 0;
  int iz, iy, ix;
  int offset[3];
  int dep, a, diff;
  int gldm_idx;
  for (iz = 0; iz < size[0]; iz ++) // Iterate over all voxels in a row - column - slice order
  {
    for (iy = 0; iy < size[1]; iy ++)
    {
      for (ix = 0; ix < size[2]; ix ++)
      {
        if (mask[i])  // Check if the current voxel is part of the segmentation
        {
          dep = 0;
          for (a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
          {
            offset[0] = angles[a * 3];
            offset[1] = angles[a * 3 + 1];
            offset[2] = angles[a * 3 + 2];

            // Check whether the neighbour index is not out of range (i.e. part of the image)
            if (iz + offset[0] >= 0 && iz + offset[0] < size[0] &&
                iy + offset[1] >= 0 && iy + offset[1] < size[1] &&
                ix + offset[2] >= 0 && ix + offset[2] < size[2])
            {
              j = i + offset[0] * strides[0] +
                      offset[1] * strides[1] +
                      offset[2] * strides[2];
              if (mask[j])  // Check whether neighbour voxel is part of the segmentation
              {
                diff = image[i] - image[j];
                if (diff < 0) diff *= -1;  // Get absolute difference
                if (diff <= alpha) dep++;
              }
            }
          }
          gldm_idx = dep + (image[i]-1) * (Na * 2 + 1);  // Row_idx (dep) + Col_idx (Gray level * Max dependency)
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

int get_angle_count(int *size, int *distances, int n_dim, int n_dist, char bidirectional, int force2Ddim)
{
  int n_a;
  int Na_d, Na_dd;  // Angles per distance
  int dist_idx, dim_idx;

  // First, determine the maximum distance and the number of angles to compute
  // Number of angles to compute for distance Na_d = (2d + 1)^n_dim - (2d - 1)^n_dim
  // The first term is temporarily stored in Na_d, the second in Na_dd
  n_a = 0;
  for (dist_idx = 0; dist_idx < n_dist; dist_idx++)
  {
    if (distances[dist_idx] < 1) return 0;  // invalid distance encountered

    Na_d = 1;
    Na_dd = 1;
    for (dim_idx = 0; dim_idx < n_dim; dim_idx++)
    {
        // Do not generate angles that move in the out-of-plane dimension
        if (dim_idx == force2Ddim) continue;

        // Check if the distance is within the image size for this dimension
        if (distances[dist_idx] < size[dim_idx])
        {
          // Full range possible, in this dimension, so multiply by (2d + 1) and (2d -1)
          Na_d *= (2 * distances[dist_idx] + 1);
          Na_dd *= (2 * distances[dist_idx] - 1);
        }
        else
        {
          // Limited range possible, so multiply by (2 * (size - 1) + 1) (i.e. max possible distance for this size)
          Na_d *= (2 *(size[dim_idx] - 1) + 1);
          Na_dd *= (2 *(size[dim_idx] - 1) + 1);
        }
    }
    n_a += (Na_d - Na_dd);  // Add the number of angles to be generated for this distance to the grand total
  }

  // if only single direction is needed, divide Na by 2 (if bidirectional, Na will be even, and angles is a mirrored)
  if (!bidirectional) n_a /= 2;

  return n_a;
}

int build_angles(int *size, int *distances, int n_dim, int n_dist, int force2Ddim, int n_a, int *angles)
{
  int *offset_stride;
  int max_distance, n_offsets, offset, a_dist;
  int dist_idx, dim_idx, a_idx, new_a_idx;

  max_distance = 0;  // Maximum offset specified, needed later on to generate the range of offsets
  for (dist_idx = 0; dist_idx < n_dist; dist_idx++)
  {
      if (distances[dist_idx] < 1) return 1;  // invalid distance encountered

      // Store maximum distance specified
      if (max_distance < distances[dist_idx])
        max_distance = distances[dist_idx];
  }

  n_offsets = 2 * max_distance + 1;  // e.g. for max distance = 2, offsets = {-2, -1, 0, 1, 2}; ||offsets|| = 5

  // offset_stride is used to generate the unique combinations of offsets for the angles
  // For the last dimension the stride is 1, i.e. a different offset is used for each subsequent angle
  // For the next-to-last dimension the stide is n_offsets, i.e. for each cycle through possible offsets in the last
  // dimension, the offset in this dimension changes by 1.
  // For subsequent dimensions, the stride is the previous stride multiplied by n_offsets, allowing the previous
  // dimension to cycle through all possible offsets before advancing the offset in this dimension. e.g.:
  // stride {9, 3, 1}. This reversed order ensures compatibility with Python generated angles
  // angle   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  // dim 2  -1  0  1 -1  0  1 -1  0  1 -1  0  1 -1  0  1 -1  0  1 -1  0  1 -1  0  1 -1  0  1
  // dim 1  -1 -1 -1  0  0  0  1  1  1 -1 -1 -1  0  0  0  1  1  1 -1 -1 -1  0  0  0  1  1  1
  // dim 0  -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1
  offset_stride = (int *)calloc(n_dim, sizeof(int));
  offset_stride[n_dim - 1] = 1;
  for (dim_idx = n_dim - 2; dim_idx >= 0; dim_idx--)
  {
    offset_stride[dim_idx] = offset_stride[dim_idx + 1] * n_offsets;
  }

  new_a_idx = 0;  // index used to generate new angle offset, increases during every loop
  a_idx = 0;  // index in angles array of current angle being generated, increases only when a valid angle has been generated
  while (a_idx < n_a)
  {
    a_dist = 0;  // maximum offset of the angle, corresponds to the distance this angle belongs to (infinity norm)
    // generate new angle
    for (dim_idx = 0; dim_idx < n_dim; dim_idx++)
    {
      offset = max_distance - (new_a_idx / offset_stride[dim_idx]) % n_offsets;  // {max_d, ... , -max_d}, step 1
      if ((dim_idx == force2Ddim && offset != 0) ||  // Moves in an invalid direction (out-of-plane dimension)
          offset >= size[dim_idx] ||  // Offset (positive) larger than size
          offset <= -size[dim_idx])   // Offset (negative) smaller than negative size
      {
        a_dist = -1;  // invalid angle
        break;  // no need to generate offsets for other dimensions, angle is invalid
      }
      angles[a_idx * n_dim + dim_idx] = offset;

      if (a_dist < offset) a_dist = offset;  // offset positive
      else if (a_dist < -offset) a_dist = -offset;  // offset negative
    }
    new_a_idx++;  // always advance new_a_idx, this controls the combination of offsets in generating the angle

    if (a_dist < 1) continue; // Angle is invalid, i.e a_dist = -1 (failed check) or a_dist = 0 (angle (0, 0, 0))

    // Check if the distance this angle is generated for is requested (i.e. present in distances)
    for (dist_idx = 0; dist_idx < n_dist; dist_idx++)
    {
      if (a_dist == distances[dist_idx])
      {
        a_idx++; // Angle valid, store it and move to the next
        break;  // No need to check the other distances
      }
    }
  }

  free(offset_stride);
  return 0;
}
