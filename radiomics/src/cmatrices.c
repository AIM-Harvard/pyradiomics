#include <stdlib.h>
#include "cmatrices.h"

int calculate_glcm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, double *glcm, int Ng)
{
  /* Calculate GLCM: Count the number of voxels with gray level i is neighboured by a voxel with gray level j in
  *  direction and distance specified by a. Returns an asymmetrical GLCM matrix for each angle/distance
  *  defined in angles.
  */

  // Index and size variables of the image
  size_t Ni;  // Size of the entire image array
  size_t i, j;  // Iterator variables (image)
  size_t* cur_idx = (size_t *)malloc(sizeof *cur_idx * Nd);  // Temporary array to store current index by dimension

  size_t a, d;  // Iterator variables (angles, dimensions)

  // Output matrix variables
  size_t glcm_idx, glcm_idx_max = Ng * Ng * Na;  // Index and max index of the texture array

  // Calculate size of image array, and set i at lower bound of bounding box
  Ni = size[0];
  i = bb[0] * strides[0];
  for (d = 1; d < Nd; d++)
  {
    i += bb[d] * strides[d];
    Ni *= size[d];
  }

  // Loop over all voxels in the image
  for ( ; i < Ni; i++)
  {
    // Calculate the current index in each dimension
    for (d = Nd - 1; d > 0; d--)  // Iterate in reverse direction to handle strides from small to large
    {
      cur_idx[d] = (i % strides[d - 1]) / strides[d];
      if (cur_idx[d] > bb[Nd + d])
      {
        // Set the i to the lower bound of the bounding box
        // size[d] - cur_idx[d] ensures an overflow, causing the index in current dimension to be 0
        // Then, add bb[d] to ensure it is set to the lower bound of the bounding box
        i += (size[d] - cur_idx[d] + bb[d]) * strides[d];
        cur_idx[d] = bb[d];  // Set cur_idx[d] to reflect the change to i
      }
      else if (cur_idx[d] < bb[d])
      {
        i += (bb[d] - cur_idx[d]) * strides[d];
        cur_idx[d] = bb[d];  // Set cur_idx[d] to reflect the change to i
      }
    }

    cur_idx[0] = i / strides[0];
    if (cur_idx[0] > bb[Nd])  // No need to check < bb[d], as initialization sets this at bb[d]
      break; // Out-of-range in first dimension: end of bounding box reached

    if (mask[i])
    {
      // Loop over all angles to get the neighbours
      for (a = 0; a < Na; a++)
      {
        j = i;  // Start at current center voxel
        for (d = 0; d < Nd; d++)
        {
          // Check if the current offset does not go out-of-range
          if (cur_idx[d] + angles[a * Nd + d] < bb[d] ||  cur_idx[d] + angles[a * Nd + d] > bb[Nd + d])
          {
            // Set to i to signal out-of-range below
            // (normally j would never be equal to i, as an angle has at least 1 non-zero offset)
            j = i;
            break;
          }
          j += angles[a * Nd + d] * strides[d];
        }

        // If the neighbor voxel is not out of range (signalled by setting j=i) and part of the ROI (mask[j]),
        // increment the corresponding element in the GLCM
        if (j != i && mask[j])
        {
          glcm_idx = a + (image[j]-1) * Na + (image[i]-1) * Na * Ng;
          if (image[i] <= 0 || image[j] <= 0 || glcm_idx >= glcm_idx_max)
          {
            free(cur_idx);
            return 0; // Index out of range
          }
          glcm[glcm_idx] ++;
        }
      }
    }
  }
  free(cur_idx);
  return 1;
}

int calculate_glszm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, int *tempData, int Ng, int Ns, int Nvox)
{
  /* Calculate the GLSZM: Count the number of connected zones with gray level i and size j in the ROI. Uses the angles
   * to find neighbours and pushes neighbours with the same gray level of the current zone onto a stack. Next, the last
   * entry to the stack is popped and neighbours are again evaluated and pushed. If no the stack is empty, the region is
   * complete. The number of times a voxel has been pushed / popped from the stack is the size of the zone. Voxels that
   * have been assigned to a zone are marked to prevent re-processing.
   *
   * This function calculates the zones and stores them in a temporary output array, combining the size with the gray
   * level. A separate function the fills the GLSZM array from this tempData array. This allows memory conservation,
   * as the GLSZM only has to be instantiated to allow the found sizes, rather than the potential sizes based on voxel
   * count alone.
   */

  // Index and size variables of the image
  size_t Ni;  // Size of the entire image array
  size_t i, j, k;  // Iterator variables (image)
  size_t* cur_idx = (size_t *)malloc(sizeof *cur_idx * Nd);  // Temporary array to store current index by dimension

  // Stack to hold indices of a growing region
  size_t *regionStack;
  size_t stackTop = 0;
  size_t *processedStack = NULL;
  size_t processed_idx = 0;

  size_t a, d;  // Iterator variables (angles, dimensions)

  // Output matrix variables
  int gl, region;
  int maxSize = 0;
  size_t temp_idx = 0;
  size_t temp_idx_max = Ns * 2;

  regionStack = (size_t *)malloc(sizeof *regionStack * Ns);

  // If processing multiple voxels, use a processedStack to keep track of processed voxels.
  // These need to be reset after processing to allow for reprocessing in the next kernel(s)
  if (Nvox > 1)
    processedStack = (size_t *)malloc(sizeof *processedStack * Ns);

  // Calculate size of image array, and set i at lower bound of bounding box
  Ni = size[0];
  i = bb[0] * strides[0];
  for (d = 1; d < Nd; d++)
  {
    i += bb[d] * strides[d];
    Ni *= size[d];
  }

  // Loop over all voxels in the image
  for ( ; i < Ni; i++)
  {
    // Calculate the current index in each dimension (except dimension 0, handled below)
    for (d = Nd - 1; d > 0; d--)  // Iterate in reverse direction to handle strides from small to large
    {
      cur_idx[d] = (i % strides[d - 1]) / strides[d];
      if (cur_idx[d] > bb[Nd + d])
        // Set the i to the lower bound of the bounding box
        // size[d] - cur_idx[d] ensures an overflow, causing the index in current dimension to be 0
        // Then, add bb[d] to ensure it is set to the lower bound of the bounding box
        i += (size[d] - cur_idx[d] + bb[d]) * strides[d];
      else if (cur_idx[d] < bb[d])
        i += (bb[d] - cur_idx[d]) * strides[d];
    }

    if (i / strides[0] > bb[Nd])  // No need to check < bb[d], as initialization sets this at bb[d]
      break; // Out-of-range in first dimension: end of bounding box reached

    // Check if the current voxel is part of the segmentation and unprocessed
    if (mask[i])
    {
      // Store current gray level, needed for region growing and determines index in GLSZM.
      gl = image[i];

      // Instantiate variable to hold region size at 0. Region increases for every found voxel.
      region = 0;

      // Voxel-based: Add the current voxel to the processed stack to reset later.
      if (processedStack)
      {
        if (processed_idx >= Ns)  // index out of range
        {
          free(cur_idx);
          free(regionStack);
          return -1;
        }
        // ++ after processed_idx --> use first, then increase
        // (i.e. returned value = x, processed_idx value after statement = x + 1)
        processedStack[processed_idx++] = i;
      }

      // Start growing the region
      regionStack[stackTop++] = i; // Add the current voxel to the stack as 'starting point'
      mask[i] = 0;  // Mark current voxel as 'processed'

      while (stackTop > 0)
      {
        k = regionStack[--stackTop];  // Get the next voxel to process, on first iteration, this equals i

        // Increment region size, as number of loops corresponds to number of voxels in current region
        region++;

        // Calculate the current index in each dimension. No checks with bb are needed, as k is either
        // equal to i (subject to bb checks) or to j (subject to bb checks)
        cur_idx[0] = k / strides[0];
        for (d = 1; d < Nd; d++)
          cur_idx[d] = (k % strides[d - 1]) / strides[d];

        // Generate neighbours for current voxel
        for (a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
        {
          j = k;  // Start at current center voxel
          for (d = 0; d < Nd; d++)
          {
            // Check if the current offset does not go out-of-range
            if (cur_idx[d] + angles[a * Nd + d] < bb[d] ||  cur_idx[d] + angles[a * Nd + d] > bb[Nd + d])
            {
              // Set to k to signal out-of-range below
              // (normally j would never be equal to k, as an angle has at least 1 non-zero offset)
              j = k;
              break;
            }
            j += angles[a * Nd + d] * strides[d];
          }

          // If the neighbor voxel is not out of range (> 0) and part of the ROI (mask[j]),
          // increment the corresponding element in the GLCM
          if (j != k && mask[j] && (image[j] == gl))
          {
            // Voxel-based: Add the current voxel to the processed stack to reset later.
            if (processedStack)
            {
              if (processed_idx >= Ns)  // index out of range
              {
                free(cur_idx);
                free(regionStack);
                return -1;
              }
              // ++ after processed_idx --> use first, then increase
              // (i.e. returned value = x, processed_idx value after statement = x + 1)
              processedStack[processed_idx++] = j;
            }

            // Push the voxel index to the stack for further processing
            regionStack[stackTop++] = j;
            // Voxel belongs to current region, mark it as 'processed'
            mask[j] = 0;
          }
        } // next a
      }  // while (stackTop > -1)

      if (temp_idx >= temp_idx_max)  // index out of range
      {
        free(cur_idx);
        free(regionStack);
        return -1;
      }
      // Keep track of the largest region encountered, used to instantiate the GLSZM matrix later
      if (region > maxSize) maxSize = region;

      // Store the region size and gray level in the temporary output matrix
      tempData[(temp_idx * 2)] = gl;
      tempData[((temp_idx * 2) + 1)] = region;

      temp_idx++;
    }
  }
  free(cur_idx);
  free(regionStack);

  // Reset all processed voxels (needed when computing voxel-based)
  if (processedStack)
  {
    // -- before processed_idx --> then increase first, then use
    // (i.e. returned value = x - 1, processed_idx value after statement = x - 1)
    while (processed_idx > 0)
      mask[processedStack[--processed_idx]] = 1;
    free(processedStack);
  }

  if (temp_idx >= temp_idx_max)
    return -1; // index out of range
  tempData[(temp_idx * 2)] = -1; // Set the first element after last region to -1 to stop the loop in fill_glszm

  return maxSize;
}

int fill_glszm(int *tempData, double *glszm, int Ng, int maxRegion)
{
  /* This function fills the GLSZM using the zones described in the tempData. See calculate_glszm() for more details.
   */
  size_t i = 0;
  size_t glszm_idx, glszm_idx_max = Ng * maxRegion;  // Index and max index of the texture array

  while(tempData[i * 2] > -1)
  {
    glszm_idx = (tempData[i * 2] - 1) * maxRegion + tempData[i * 2 + 1] - 1;
    if (tempData[i * 2] <= 0 || glszm_idx >= glszm_idx_max) return 0; // Index out of range

    glszm[glszm_idx]++;
    i++;
  }
  return 1;
}

int calculate_glrlm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, double *glrlm, int Ng, int Nr)
{
  /* Calculate the GLRLM: Count the number of runs (=consecutive voxels with the same gray level) with gray level i and
   * run length j along angle a.
   *
   * The GLRLM is calculated by angle: First, it detects in which dimensions the current angle has a non-zero offset,
   * and what direction is the start (if the angle is positive, the start is at index 0, if negative, at size - 1).
   * Then the function loops over all the voxels in the image and checks if the voxel represents a starting position for
   * any of the moving dimensions. If not, the index is skipped forward to a point that is a starting position.
   * Finally, the run is instantiated by starting at the current index and advancing until the index goes out-of-range
   * in a dimension. During the run, the code tracks the number of consecutive voxels the define a run and updates the
   * GLRLM accordingly.
   */
  // Index and size variables of the image
  size_t Ni, start_i;  // Size of the entire image array
  size_t i, j;  // Iterator variables (image)
  size_t cur_idx;  // Only need a single int for index, as each the index is calculated and check separately for each dimension
  int start_voxel_valid;  // boolean indicating whether the current voxel is a valid starting voxel

  size_t a, d, md;  // Iterator variables (angles, dimensions)
  size_t cnt_mDim;  // Variable to hold the number of moving dims
  int multiElement;  // Variable to check whether the current angle only yields runs of length 1

  // Variables to track the non-zero offsets of the current angle (=moving dimensions)
  // and to help define the start voxels for the runs
  size_t* mDims = (size_t *)malloc(sizeof *mDims * Nd);  // Array to hold mapping to moving dimensions
  size_t* mDim_start = (size_t *)malloc(sizeof *mDim_start * Nd); // Array to hold start positions for moving dimensions (0 or size - 1)

  // Output matrix variables
  int gl, rl, elements;
  size_t glrlm_idx, glrlm_idx_max = Ng * Nr * Na;  // Index and max index of the texture array

  // Calculate size of image array, and calculate index of lower bound of bounding box (`start_i`)
  Ni = size[0];
  start_i = bb[0] * strides[0];
  for (d = 1; d < Nd; d++)
  {
    start_i += bb[d] * strides[d];
    Ni *= size[d];
  }

  for (a = 0; a < Na; a++)  // Iterate over angles to get the neighbours
  {
    multiElement = 0;
    // First lookup and count the number of dimensions where the angle != 0 (i.e. "Moving dimensions")
    // Moreover, check if we need to start at 0 (angle > 0) or at the end (size[d] - 1, angle < 0)
    cnt_mDim = 0;
    for (d = 0; d < Nd; d++)
    {
      if (angles[a * Nd + d] != 0)
      {
        if (angles[a * Nd + d] > 0)
          mDim_start[cnt_mDim] = bb[d];
        else
          mDim_start[cnt_mDim] = bb[Nd + d];
        mDims[cnt_mDim] = d;
        cnt_mDim++;
      }
    }

    // Then, iterate over the image (with the goal of getting all "start positions", i.e. from where we start the run)
    for (i = start_i; i < Ni; i++)
    {
      // Calculate the current index in each dimension and ensure it is within the bounding box.
      // (except dimension 0, handled below)
      for (d = Nd - 1; d > 0; d--)  // Iterate in reverse direction to handle strides from small to large
      {
        cur_idx = (i % strides[d - 1]) / strides[d];
        if (cur_idx > bb[Nd + d])
          // Set the i to the lower bound of the bounding box
          // size[d] - cur_idx[d] (i.e. add the difference between current index and size to i  --> advance along this
          // dimension to index 'size'. As the max index in each dimension is equal to size - 1,
          // this ensures an overflow, causing the index in current dimension to be 0 and in a lower dimension to be +1.
          // Then, add bb[d] to ensure it is set to the lower bound of the bounding box in the current dimension
          i += (size[d] - cur_idx + bb[d]) * strides[d];
        else if (cur_idx < bb[d])
          i += (bb[d] - cur_idx) * strides[d];
      }

      cur_idx = i / strides[0];
      if (i / strides[0] > bb[Nd])  // No need to check < bb[d], as initialization sets this at bb[d]
        break; // Out-of-range in first dimension: end of bounding box reached

      // loop over moving dimensions to check if the current voxel is a valid starting voxel
      // a voxel is a starting voxel if it's index for a specific moving dimension matches the start index for
      // that dimension. If this is the case for ANY of the moving dimensions, the voxel is indeed a starting voxel
      start_voxel_valid = 0;
      for (md = 0; md < cnt_mDim; md++)
      {
        d = mDims[md];
        // Get the current index for the current moving dimension mDims[d]
        if (d == 0)
          cur_idx = i / strides[d];
        else
          cur_idx = (i % strides[d - 1]) / strides[d];

        if (cur_idx == mDim_start[md])
        {
          start_voxel_valid = 1;
          break;  // No need to check the rest, the voxel is valid
        }
      }

      if (!start_voxel_valid)
      {
        // Oh oh, current voxel is not a starting position for any of the moving dimensions!
        // Skip to a valid index by ensuring the fastest changing moving dimension is set to a valid start position.
        md = cnt_mDim - 1;
        d = mDims[md]; // Get the last moving dimension (i.e. the moving dimension with the smallest stride)

        // Advance i in the last moving dimension to a valid start position. Do this by calculating the difference
        // between the current position (cur_idx, as the loop above ends with the last moving dimension) and the
        // intended start position.
        // Take the modulus with the size to ensure this change is always forward.
        // Add size[d] to ensure the operation returns the modulus, not the remainder (% operator)
        i += ((mDim_start[md] - cur_idx + size[d]) % size[d]) * strides[d]; // Skip the rest of the voxels in this moving dimension

        if (d > 1)
        {
          // Update all lower dimensions if necessary (ensuring it all stays within the bounding box)
          d--;  // Don't handle the current dimension (=last moving dimension), we just did that above...
          for ( ; d > 0 ; d--)  // don't handle d = 0, that is done below...
          {
            cur_idx = (i % strides[d - 1]) / strides[d];
            if (cur_idx > bb[Nd + d])
            {
              // Set the i to the lower bound of the bounding box
              // size[d] - cur_idx[d] ensures an overflow, causing the index in current dimension to be 0
              // Then, add bb[d] to ensure it is set to the lower bound of the bounding box
              i += (size[d] - cur_idx + bb[d]) * strides[d];
            }
            else if (cur_idx < bb[d])
              i += (bb[d] - cur_idx) * strides[d];
          }
        }

        if (i / strides[0] > bb[Nd])
          break; // Out-of-range in first dimension: end of bounding box reached
      }

      // Run Forest, Run! Start at the current index and advance using the angle until exiting the image.
      j = i;
      gl = -1;
      rl = 0;
      elements = 0;

      // Run the body at least once, without checking if j==i.
      // On the first execute j==i, but afterwards is changed. If j==i at the end of this do...while block, it is
      // to signal that j has gone out-of-range and the run is done.
      do
      {
        // Check if j belongs to the mask (i.e. part of some run)
        if (mask[j])
        {
          elements++;  // count the number of voxels in this run
          if (gl == -1)  // No run initialized, start a new one
            gl = image[j];
          else if (image[j] == gl)  // j is part of the current run, increase the length
            rl++;
          else  // j is not part of the run, end current run and start new one
          {
            glrlm_idx = a + rl * Na + (gl - 1) * Na * Nr;
            if (gl <= 0 || glrlm_idx >= glrlm_idx_max)
            {
              free(mDims);
              free(mDim_start);
              return 0;
            }
            glrlm[glrlm_idx]++;

            gl = image[j];
            rl = 0;
          }
        }
        else if (gl > -1)  // end current run
        {
          glrlm_idx = a + rl * Na + (gl - 1) * Na * Nr;
          if (gl <= 0 || glrlm_idx >= glrlm_idx_max)
          {
            free(mDims);
            free(mDim_start);
            return 0;
          }
          glrlm[glrlm_idx]++;

          gl = -1;
          rl = 0;
        }

        // Advance to the next voxel:
        // Set j to the next voxel using the current angle. If j goes out-of-range, it is set to -1,
        // causing the while(j > -1) to exit (end of the runs for the current index).
        for (md = 0; md < cnt_mDim; md++)
        {
          d = mDims[md];
          // Calculate the current index in each dimension
          if (d == 0)
            cur_idx = j / strides[d];
          else
            cur_idx = (j % strides[d - 1]) / strides[d];

          // Check if that is within the range [0, size[d])
          if (cur_idx + angles[a * Nd + d] < bb[d] ||  cur_idx + angles[a * Nd + d] > bb[Nd + d])
          {
            j = i;  // Set j to i to signal the while loop to exit too.
            break;  // Out of range! Run is done, so exit the loop
          }

          j += angles[a * Nd + d] * strides[d];
        }  // next md
      }
      while (j != i);

      if (gl > -1)  // end current run (this is the case when the last voxel of the run was included in the mask)
      {
        glrlm_idx = (gl - 1) * Na * Nr + rl * Na + a;
        if (gl <= 0 || glrlm_idx >= glrlm_idx_max)
        {
          free(mDims);
          free(mDim_start);
          return 0;
        }
        glrlm[glrlm_idx]++;
      }

      if (elements > 1) // multiple elements found in this run. Ensure multiElement is set to true;
        multiElement = 1;
    }  // next i

    if (!multiElement) // Segmentation is 2D for this angle, remove it.
    {
      // set elements at runlength index 0 to 0 for all gray levels in this angle, other runlength indices are
      // already 0.
      for (glrlm_idx = 0; glrlm_idx < Ng; glrlm_idx++)
        glrlm[glrlm_idx * Nr * Na + a] = 0;
    }
  }  // next a

  free(mDims);
  free(mDim_start);

  return 1;
}

int calculate_ngtdm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, double *ngtdm, int Ng)
{
  /* Calculate the NGTDM: For each voxel, calculate the absolute difference between the center voxel and the average of
   * its neighbours. Then, add this difference to a grand total (for the gray level i of the center voxel)
   *
   */
  // Index and size variables of the image
  size_t Ni;  // Size of the entire image array
  size_t i, j;  // Iterator variables (image)
  size_t* cur_idx = (size_t *)malloc(sizeof *cur_idx * Nd);  // Temporary array to store current index by dimension

  size_t a, d;  // Iterator variables (angles, dimensions)

  // Output matrix variables
  int gl;
  double count, sum, diff;
  size_t ngtdm_idx, ngtdm_idx_max = Ng * 3;  // Index and max index of the texture array

  // Fill gray levels (empty slices gray levels are later deleted in python)
  for (gl = 0; gl < Ng; gl++)
  {
      ngtdm[gl*3 + 2] = gl + 1;
  }
  /* Calculate matrix: for each gray level, element 0 describes the number of voxels with gray level i and
   * element 1 describes the sum of all differences between voxels with gray level i and their neighbourhood.
   * element 2 is set above and contains the corresponding gray level (gray levels not present in the ROI are removed
   * later on).
   */

  // Calculate size of image array, and set i at lower bound of bounding box
  Ni = size[0];
  i = bb[0] * strides[0];
  for (d = 1; d < Nd; d++)
  {
    i += bb[d] * strides[d];
    Ni *= size[d];
  }

  // Loop over all voxels in the image
  for ( ; i < Ni; i++)
  {
    // Calculate the current index in each dimension (except dimension 0, handled below)
    for (d = Nd - 1; d > 0; d--)  // Iterate in reverse direction to handle strides from small to large
    {
      cur_idx[d] = (i % strides[d - 1]) / strides[d];
      if (cur_idx[d] > bb[Nd + d])
      {
        // Set the i to the lower bound of the bounding box
        // size[d] - cur_idx[d] ensures an overflow, causing the index in current dimension to be 0
        // Then, add bb[d] to ensure it is set to the lower bound of the bounding box
        i += (size[d] - cur_idx[d] + bb[d]) * strides[d];
        cur_idx[d] = bb[d];  // Set cur_idx[d] to reflect the change to i
      }
      else if (cur_idx[d] < bb[d])
      {
        i += (bb[d] - cur_idx[d]) * strides[d];
        cur_idx[d] = bb[d];  // Set cur_idx[d] to reflect the change to i
      }
    }

    cur_idx[0] = i / strides[0];
    if (cur_idx[0] > bb[Nd])  // No need to check < bb[d], as initialization sets this at bb[d]
      break; // Out-of-range in first dimension: end of bounding box reached

    if (mask[i])
    {
      count = 0;
      sum = 0;

      // Loop over all angles to get the neighbours
      for (a = 0; a < Na; a++)
      {
        j = i;  // Start at current center voxel
        for (d = 0; d < Nd; d++)
        {
          // Check if the current offset does not go out-of-range
          if (cur_idx[d] + angles[a * Nd + d] < bb[d] ||  cur_idx[d] + angles[a * Nd + d] > bb[Nd + d])
          {
            // Set to i to signal out-of-range below
            // (normally j would never be equal to i, as an angle has at least 1 non-zero offset)
            j = i;
            break;
          }
          j += angles[a * Nd + d] * strides[d];
        }

        // If the neighbor voxel is not out of range (signalled by setting j=i) and part of the ROI (mask[j]),
        // increment the corresponding element in the GLCM
        if (j != i && mask[j])
        {
           count++;
           sum += image[j];
        }
      }
      if (count == 0)
        diff = 0;
      else
        diff = (double)image[i] - (sum  / count);

      if (diff < 0)
        diff *= -1;  // Get absolute difference

      ngtdm_idx = (image[i]-1) * 3;
      if (image[i] <= 0 || ngtdm_idx >= ngtdm_idx_max)
      {
        free(cur_idx);
        return 0; // Index out of range
      }
      ngtdm[ngtdm_idx]++;  // Increase the count for this gray level
      ngtdm[ngtdm_idx + 1] += diff;  // Add the absolute difference to the total for this gray level
    }
  }

  free(cur_idx);
  return 1;
}

int calculate_gldm(int *image, char *mask, int *size, int *bb, int *strides, int *angles, int Na, int Nd, double *gldm, int Ng, int alpha)
{
  /* Calculate GLDM: Count the number of voxels with gray level i, that have j dependent neighbours.
  *  A voxel is considered dependent if the absolute difference between the center voxel and the neighbour <= alpha
  */

  // Index and size variables of the image
  size_t Ni;  // Size of the entire image array
  size_t i, j;  // Iterator variables (image)
  size_t* cur_idx = (size_t *)malloc(sizeof *cur_idx * Nd);  // Temporary array to store current index by dimension

  size_t a, d;  // Iterator variables (angles, dimensions)

  // Output matrix variables
  int dep, diff;
  size_t gldm_idx, gldm_idx_max = Ng * (Na * 2 + 1);  // Index and max index of the texture array

  // Calculate size of image array, and set i at lower bound of bounding box
  Ni = size[0];
  i = bb[0] * strides[0];
  for (d = 1; d < Nd; d++)
  {
    i += bb[d] * strides[d];
    Ni *= size[d];
  }

  // Loop over all voxels in the image
  for ( ; i < Ni; i++)
  {
    // Calculate the current index in each dimension (except dimension 0, handled below)
    for (d = Nd - 1; d > 0; d--)  // Iterate in reverse direction to handle strides from small to large
    {
      cur_idx[d] = (i % strides[d - 1]) / strides[d];
      if (cur_idx[d] > bb[Nd + d])
      {
        // Set the i to the lower bound of the bounding box
        // size[d] - cur_idx[d] ensures an overflow, causing the index in current dimension to be 0
        // Then, add bb[d] to ensure it is set to the lower bound of the bounding box
        i += (size[d] - cur_idx[d] + bb[d]) * strides[d];
        cur_idx[d] = bb[d];  // Set cur_idx[d] to reflect the change to i
      }
      else if (cur_idx[d] < bb[d])
      {
        i += (bb[d] - cur_idx[d]) * strides[d];
        cur_idx[d] = bb[d];  // Set cur_idx[d] to reflect the change to i
      }
    }

    cur_idx[0] = i / strides[0];
    if (cur_idx[0] > bb[Nd])  // No need to check < bb[d], as initialization sets this at bb[d]
      break; // Out-of-range in first dimension: end of bounding box reached

    if (mask[i])
    {
      dep = 0;

      // Loop over all angles to get the neighbours
      for (a = 0; a < Na; a++)
      {
        j = i;  // Start at current center voxel
        for (d = 0; d < Nd; d++)
        {
          // Check if the current offset does not go out-of-range
          if (cur_idx[d] + angles[a * Nd + d] < bb[d] ||  cur_idx[d] + angles[a * Nd + d] > bb[Nd + d])
          {
            // Set to i to signal out-of-range below
            // (normally j would never be equal to i, as an angle has at least 1 non-zero offset)
            j = i;
            break;
          }
          j += angles[a * Nd + d] * strides[d];
        }

        // If the neighbor voxel is not out of range (signalled by setting j=i) and part of the ROI (mask[j]),
        // increment the corresponding element in the GLCM
        if (j != i && mask[j])
        {
          diff = image[i] - image[j];
          if (diff < 0) diff *= -1;  // Get absolute difference
          if (diff <= alpha) dep++;
        }
      }
      gldm_idx = dep + (image[i]-1) * (Na * 2 + 1);  // Row_idx (dep) + Col_idx (Gray level * Max dependency)
      if (image[i] <= 0 || gldm_idx >= gldm_idx_max)
      {
        free(cur_idx);
        return 0; // Index out of range
      }
      gldm[gldm_idx] ++;
    }
  }

  free(cur_idx);
  return 1;
}

int get_angle_count(int *size, int *distances, int Nd, int Ndist, char bidirectional, int force2Ddim)
{
  /* Calculate the number of angles that need to be generated for the requested distances, taking into account the size
   * of the image and whether a dimension has been excluded (=force2Ddim).
   *
   * First, determine the maximum distance and the number of angles to compute
   * Number of angles to compute for distance Na_d = (2d + 1)^Nd - (2d - 1)^Nd
   * The first term is temporarily stored in Na_d, the second in Na_dd
   */
  int Na;  // Grand total of angles to generate
  int Na_d, Na_dd;  // Angles per distance Na_d = N angles for distance [0, n] and Na_dd = N angles for distance [0, n -1]
  int dist_idx, dim_idx;

  Na = 0;
  for (dist_idx = 0; dist_idx < Ndist; dist_idx++)
  {
    if (distances[dist_idx] < 1) return 0;  // invalid distance encountered

    Na_d = 1;
    Na_dd = 1;
    for (dim_idx = 0; dim_idx < Nd; dim_idx++)
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
          // This is the same for both Na_d and Na_dd, as the upper bound is the size, not the distance.
          // The multiplication is still needed, as Na_d and Na_dd diverge more strongly, and are subtracted later on.
          Na_d *= (2 *(size[dim_idx] - 1) + 1);
          Na_dd *= (2 *(size[dim_idx] - 1) + 1);
        }
    }
    Na += (Na_d - Na_dd);  // Add the number of angles to be generated for this distance to the grand total
  }

  // if only single direction is needed, divide Na by 2
  // (if bidirectional, Na will be even, and angles is a mirrored array)
  if (!bidirectional) Na /= 2;

  return Na;
}

int build_angles(int *size, int *distances, int Nd, int Ndist, int force2Ddim, int Na, int *angles)
{
  /* Generate the angles. One index, new_a_idx controls the combination of offsets, and is always increased.
   * Each new generated angle is checked if it is valid (i.e. at least 1 dimension non-zero, no offset larger than
   * the size and, if specified, only 0 offset in the force2Ddim). During generation, the offsets for the angle are
   * stored in the angles array, but a_idx is only increased if the angle is valid and it's maximum offset (defining
   * the infinity norm distance the angle belongs to) is specified in the distances array. That way, invalid angles are
   * overwritten by the next generated angle. This loop continues until the specified number of angles are generated.
   *
   * N.B. This function should only be used with the results obtained from `get_angle_count`, and values for `size`,
   * `distances`, `Nd`, `Ndist` and `force2Ddim` should be identical to those used in `get_angle_count`.
   */

  int *offset_stride;
  int max_distance, n_offsets, offset, a_dist;
  int dist_idx, dim_idx, a_idx, new_a_idx;

  max_distance = 0;  // Maximum offset specified, needed later on to generate the range of offsets
  for (dist_idx = 0; dist_idx < Ndist; dist_idx++)
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
  offset_stride = (int *)malloc(sizeof *offset_stride * Nd);
  offset_stride[Nd - 1] = 1;
  for (dim_idx = Nd - 2; dim_idx >= 0; dim_idx--)
    offset_stride[dim_idx] = offset_stride[dim_idx + 1] * n_offsets;

  new_a_idx = 0;  // index used to generate new angle offset, increases during every loop
  a_idx = 0;  // index in angles array of current angle being generated, increases only when a valid angle has been generated
  while (a_idx < Na)
  {
    a_dist = 0;  // maximum offset of the angle, corresponds to the distance this angle belongs to (infinity norm)
    // generate new angle
    for (dim_idx = 0; dim_idx < Nd; dim_idx++)
    {
      offset = max_distance - (new_a_idx / offset_stride[dim_idx]) % n_offsets;  // {max_d, ... , -max_d}, step 1
      if ((dim_idx == force2Ddim && offset != 0) ||  // Moves in an invalid direction (out-of-plane dimension)
          offset >= size[dim_idx] ||  // Offset (positive) larger than size
          offset <= -size[dim_idx])   // Offset (negative) smaller than negative size
      {
        a_dist = -1;  // invalid angle
        break;  // no need to generate offsets for other dimensions, angle is invalid
      }
      angles[a_idx * Nd + dim_idx] = offset;

      if (a_dist < offset)
        a_dist = offset;  // offset positive
      else if (a_dist < -offset)
        a_dist = -offset;  // offset negative
    }
    new_a_idx++;  // always advance new_a_idx, this controls the combination of offsets in generating the angle

    if (a_dist < 1) continue; // Angle is invalid, i.e a_dist = -1 (failed check) or a_dist = 0 (angle (0, 0, 0))

    // Check if the distance this angle is generated for is requested (i.e. present in distances)
    for (dist_idx = 0; dist_idx < Ndist; dist_idx++)
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
