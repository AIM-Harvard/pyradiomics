#include "cshape.h"
#include <math.h>
#include <stdlib.h>

int *generate_angles(int *size, int *strides, int *a_strides, int *Na, int *mDim);

// Declare the look-up tables, these are filled at the bottom of this code file.
static const int gridAngles[8][3];
//static const int edgeTable[128];  // Not needed in this implementation
static const int triTable[128][16];
static const double vertList[12][3];

double calculate_surfacearea(char *mask, int *size, int *strides, double *spacing)
{
  int iz, iy, ix, i, t, d;  // iterator indices
  unsigned char cube_idx;  // cube identifier, 8 bits signifying which corners of the cube belong to the segmentation
  int a_idx;  // Angle index (8 'angles', one pointing to each corner of the marching cube
  double sum;
  double surfaceArea = 0;  // Total surface area
  double a[3], b[3], c[3];  // 2 points of the triangle, relative to the third, and the cross product vector

  // Iterate over all voxels, do not include last voxels in the three dimensions, as the cube includes voxels at pos +1
  for (iz = 0; iz < (size[0] - 1); iz++)
  {
    for (iy = 0; iy < (size[1] - 1); iy++)
    {
      for (ix = 0; ix < (size[2] - 1); ix++)
      {
        // Get current cube_idx by analyzing each point of the current cube
        cube_idx = 0;
        for (a_idx = 0; a_idx < 8; a_idx++)
        {
          i = (iz + gridAngles[a_idx][0]) * strides[0] +
              (iy + gridAngles[a_idx][1]) * strides[1] +
              (ix + gridAngles[a_idx][2]) * strides[2];

          if (mask[i]) cube_idx |= (1 << a_idx);
        }

        // Isosurface is symmetrical around the midpoint, flip the number if > 128
        // This enables look-up tables to be 1/2 the size.
        if (cube_idx & 0x80) cube_idx ^= 0xff;

        // Exlcude cubes entirely outside or inside the segmentation. Also exclude potential invalid values (>= 128).
        if (cube_idx > 0 && cube_idx < 128)
        {
          t = 0;
          while (triTable[cube_idx][t*3] >= 0) // Exit loop when no more triangles are present (element at index = -1)
          {
            for (d = 0; d < 3; d++)
            {
                a[d] = vertList[triTable[cube_idx][t*3 + 1]][d] - vertList[triTable[cube_idx][t*3]][d];
                b[d] = vertList[triTable[cube_idx][t*3 + 2]][d] - vertList[triTable[cube_idx][t*3]][d];

                // Factor in the spacing
                a[d] *= spacing[d];
                b[d] *= spacing[d];
            }

            // Compute the cross-product
            c[0] = (a[1] * b[2]) - (b[1] * a[2]);
            c[1] = (a[2] * b[0]) - (b[2] * a[0]);
            c[2] = (a[0] * b[1]) - (b[0] * a[1]);

            // Get the square
            c[0] = c[0] * c[0];
            c[1] = c[1] * c[1];
            c[2] = c[2] * c[2];

            // Compute the surface, which is equal to 1/2 magnitude of the cross product, where
            // The magnitude is obtained by calculating the euclidean distance between (0, 0, 0)
            // and the location of c
            sum = c[0] + c[1] + c[2];
            sum = sqrt(sum);
            sum = 0.5 * sum;
            surfaceArea += sum;
            t++;
          }
        }
      }
    }
  }
  return surfaceArea;
}

int calculate_diameter(char *mask, int *size, int *strides, double *spacing, int Ns, double *diameters)
{
  int iz, iy, ix, i, j;
  int a_idx, d_idx;
  int Na, mDim;
  int *angles;
  int a_strides[3];

  int *stack;
  int stack_top = -1;

  int idx, jz, jy, jx;
  double dz, dy, dx;
  double distance;

  angles = generate_angles(size, strides, a_strides, &Na, &mDim);
  stack = (int *)calloc(Ns, sizeof(int));

  // First, get all the voxels on the border
  i = 0;
  // Iterate over all voxels in a row - column - slice order
  // As the mask is padded with 0's, no voxels on the edge of the image are part of the mask, so skip those...
  for (iz = 1; iz < size[0] - 1; iz++)
  {
    for (iy = 1; iy < size[1] - 1; iy++)
    {
      for (ix = 1; ix < size[2] - 1; ix++)
      {
        i = iz * strides[0] +
            iy * strides[1] +
            ix * strides[2];
        if (mask[i])
        {
          for (a_idx = 0; a_idx < Na; a_idx++)
          {
            j = i;
            for (d_idx = 0; d_idx < mDim; d_idx++)
            {
              j += angles[a_idx * mDim + d_idx] * a_strides[d_idx];
            }

            if (mask[j] == 0)
            {
              // neighbour not part of ROI, i.e. 'i' is border voxel
              if (stack_top >= Ns) return 0;  // index out of bounds
              stack[++stack_top] = i;
              break;
            }
          }
        }
      }
    }
	}
  stack_top++; // increment by 1, so when the first item is popped, it is the last item entered

	free(angles);

  diameters[0] = 0;
  diameters[1] = 0;
  diameters[2] = 0;
  diameters[3] = 0;

  while (stack_top > 0)
  {
    // pop the last item from the stack, this prevents double processing and comparing the same voxels
    i = stack[--stack_top];
    iz = (i / strides[0]);
    iy = (i % strides[0]) / strides[1];
    ix = (i % strides[0]) % strides[1];
    for (idx = 0; idx < stack_top; idx++)  // calculate distance to all other voxels
    {
      j = stack[idx];
	    jz = (j / strides[0]);
      jy = (j % strides[0]) / strides[1];
      jx = (j % strides[0]) % strides[1];

      dz = (double)(iz - jz) * spacing[0];
      dy = (double)(iy - jy) * spacing[1];
      dx = (double)(ix - jx) * spacing[2];

      dz *= dz;
      dy *= dy;
      dx *= dx;

      distance = dz + dy + dx;
      if (iz == jz && distance > diameters[0]) diameters[0] = distance;
      if (iy == jy && distance > diameters[1]) diameters[1] = distance;
      if (ix == jx && distance > diameters[2]) diameters[2] = distance;
      if (distance > diameters[3]) diameters[3] = distance;
    }
  }
  free(stack);

  diameters[0] = sqrt(diameters[0]);
  diameters[1] = sqrt(diameters[1]);
  diameters[2] = sqrt(diameters[2]);
  diameters[3] = sqrt(diameters[3]);

  return 1;
}

int *generate_angles(int *size, int *strides, int *a_strides, int *Na, int *mDim)
{
  static int *angles;  // return value, declare static so it can be returned

  int offsets[3] = {-1, 0, 1};  // distance 1, both directions
  int a_idx, d_idx, a_offset, stride;

  // First, determine how many 'moving' dimensions there are, this determines the number of distinct angles to generate
  // Na = 3 ** NDIM(Size > 3) - 1, i.e. each dimension triples the number of angles, -1 to exclude (0, 0, 0)
  *Na = 1;
  *mDim = 0;
  for (d_idx = 0; d_idx < 3; d_idx++)  // assume mask is 3D
  {
    if (size[d_idx] > 3)  // mask is padded with 0's in all directions, i.e. bounding box size = size - 2
    {
      // Dimension is a moving dimension
      a_strides[*mDim] = strides[d_idx];
      *Na *= 3;
      (*mDim)++;
    }
  }
  (*Na)--;  // Don't generate angle for (0, 0, 0)

  // Initialize array to hold the angles
  angles = (int *)calloc(*Na * *mDim, sizeof(int));

  // Fill the angles array
  stride = 1;
  for (d_idx = 0; d_idx < *mDim; d_idx++)  // Iterate over all moving dimensions
  {
    a_offset = 0;
    for (a_idx = 0; a_idx < *Na; a_idx++)
    {
      if (a_idx == *Na / 2) a_offset = 1;  // Skip (0, 0, 0) angle

      angles[a_idx * *mDim + d_idx] = offsets[((a_idx + a_offset) / stride) % 3];
    }
    stride *= 3;  // For next dimension, multiply stride by 3 (length offsets) --> {1, 3, 9, ...}
  }
  return angles;
}

// gridAngles define the 8 corners of the marching cube, relative to the origin of the cube
static const int gridAngles[8][3] = { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 1 }, {0, 1, 0}, { 1, 0, 0 }, {1, 0, 1 }, { 1, 1, 1 }, { 1, 1, 0 } };

// edgeTable defines which edges contain intersection points, for which the exact intersection point has to be
// interpolated. However, as the intersection point is always 0.5, this can be defined beforehand, and this table is not
// needed
/*static const int edgeTable[128] = {
  0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
  0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
  0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
  0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
  0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
  0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
  0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
  0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0
};*/

// triTable defines which triangles (defined by their points as defined in vertList) are present in the cube.
// The first dimension indicates the specific cube to look up, the second dimension contains sets of 3 points (1 for
// each triangle), with the elements set to -1 after all triangles have been defined (max. no of triangles: 5)
static const int triTable[128][16] = {
  { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
  { 8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1 },
  { 3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1 },
  { 4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
  { 4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1 },
  { 9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1 },
  { 10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1 },
  { 5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
  { 5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1 },
  { 8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1 },
  { 2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
  { 2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1 },
  { 11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1 },
  { 5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1 },
  { 11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1 },
  { 11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1 },
  { 2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1 },
  { 6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
  { 3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1 },
  { 6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
  { 6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1 },
  { 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1 },
  { 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1 },
  { 3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
  { 0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1 },
  { 9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1 },
  { 8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
  { 5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1 },
  { 0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1 },
  { 6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1 },
  { 10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
  { 1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1 },
  { 0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1 },
  { 3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
  { 6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1 },
  { 9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1 },
  { 8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1 },
  { 3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1 },
  { 10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1 },
  { 10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
  { 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1 },
  { 7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
  { 2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1 },
  { 1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1 },
  { 11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1 },
  { 8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1 },
  { 0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1 },
  { 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }
};

// Vertlist represents the location of some point somewhere on an edge of the cube, relative to the origin (0, 0, 0).
// As the points on the cube are always either 0 or 1 (masked/not-masked) that other point is always halfway.
// Therefore, vertlist is constant and can be defined static (only works when the intersection point is constant
// (in this case the intersection point is always 0.5). The edge represented is defined by the gridAngle points as follows:
// { { 1, 0 }, { 2, 1 }, { 3, 2 }, { 3, 0 },
//   { 5, 4 }, { 6, 5 }, { 7, 6 }, { 7, 4 },
//   { 4, 0 }, { 5, 1 }, { 6, 2 }, { 7, 0 } }
static const double vertList[12][3] = { { 0, 0, 0.5 }, { 0, 0.5, 1 }, { 0, 1, 0.5 }, { 0, 0.5, 0 },
                                        { 1, 0, 0.5 }, { 1, 0.5, 1 }, { 1, 1, 0.5 }, { 1, 0.5, 0 },
                                        { 0.5, 0, 0 }, { 0.5, 0, 1 }, { 0.5, 1, 1 }, { 0.5, 1, 0 } };
