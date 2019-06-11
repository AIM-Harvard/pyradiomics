#include "cshape.h"
#include <math.h>
#include <stdlib.h>
#include <Python.h>


// **************************************************

// 3D Shape calculation

// **************************************************


void calculate_meshDiameter(double *vertices, size_t v_idx, double *diameters);

// Declare the look-up tables, these are filled at the bottom of this code file.
static const int gridAngles[8][3];
//static const int edgeTable[128];  // Not needed in this implementation
static const int triTable[128][16];
static const double vertList[12][3];

int calculate_coefficients(char *mask, int *size, int *strides, double *spacing,
                           double *surfaceArea, double *volume, double *diameters)
{
  int iz, iy, ix, i, t, d;  // iterator indices
  unsigned char cube_idx;  // cube identifier, 8 bits signifying which corners of the cube belong to the segmentation
  int a_idx;  // Angle index (8 'angles', one pointing to each corner of the marching cube

  static const int points_edges[2][3] = {{6, 4, 3}, {6, 7, 11}};
  size_t v_idx = 0;
  size_t v_max = 0;
  double *vertices;

  double sum;
  double a[3], b[3], c[3], ab[3];  // 3 points of the triangle, and the cross product vector
  int sign_correction;

  *surfaceArea = 0;  // Total surface area
  *volume = 0;  // Total volume

  // create a stack to hold the found vertices. For each cube, a maximum of 3 vertices are stored (with x, y and z
  // coordinates). This prevents double storing of the vertices.
  v_max = (size[0] - 1) * (size[1] - 1) * (size[2] - 1) * 9;
  vertices = (double *)calloc(v_max, sizeof(double));

  // Iterate over all voxels, do not include last voxels in the three dimensions, as the cube includes voxels at pos +1
  for (iz = 0; iz < (size[0] - 1); iz++)
  {
    for (iy = 0; iy < (size[1] - 1); iy++)
    {
      for (ix = 0; ix < (size[2] - 1); ix++)
      {
        /* Get current cube_idx by analyzing each point of the current cube
        * O - X
        * |\
        * Y Z
        *           v0
        *  p0 ------------ p1
        *   |\             |\
        *   | \ v3         | \ v1
        * v8|  \      v2   |v9\
        *   |  p3 ------------ p2
        *   |   |  v4      |   |
        *  p4 --|--------- p5  |
        *    \  |v11        \  |v10
        *  v7 \ |          v5\ |
        *      \|             \|
        *      p7 ------------ p6
        *             v6
        */
        cube_idx = 0;
        for (a_idx = 0; a_idx < 8; a_idx++)
        {
          i = (iz + gridAngles[a_idx][0]) * strides[0] +
              (iy + gridAngles[a_idx][1]) * strides[1] +
              (ix + gridAngles[a_idx][2]) * strides[2];

          if (mask[i])
            cube_idx |= (1 << a_idx);
        }

        // Isosurface is symmetrical around the midpoint, flip the number if > 128
        // This enables look-up tables to be 1/2 the size.
        // However, the sign for the volume then needs to be flipped too.
        if (cube_idx & 0x80)
        {
          cube_idx ^= 0xff;
          sign_correction = -1;
        }
        else
          sign_correction = 1;

        // ************************
        // Store vertices for diameter calculation
        // ************************

        // check if there are vertices on edges 6, 7 and 11
        // Because of the symmetry around the midpoint and the flip if cube_idx > 128, the 8th point will never appear
        // as segmented at this point. Therefore, to check if there are vertices on the adjacent edges (6, 7 and 11),
        // one only needs to check if the corresponding points (7th, 5th and 4th, respectively) are segmented.
        if (v_idx + 9 > v_max) // Overflow!
        {
          free(vertices);
          return 1;
        }

        for (t = 0; t < 3; t++)
        {
          if (cube_idx & (1 << points_edges[0][t]))
          {
            vertices[v_idx++] = (((double)iz) + vertList[points_edges[1][t]][0]) * spacing[0];
            vertices[v_idx++] = (((double)iy) + vertList[points_edges[1][t]][1]) * spacing[1];
            vertices[v_idx++] = (((double)ix) + vertList[points_edges[1][t]][2]) * spacing[2];
          }
        }

        // Exlcude cubes entirely outside or inside the segmentation (cube_idx = 0).
        if (cube_idx == 0)
          continue;

        // Process all triangles for this cube
        t = 0;
        while (triTable[cube_idx][t*3] >= 0) // Exit loop when no more triangles are present (element at index = -1)
        {
          a[0] = b[0] = c[0] = iz;
          a[1] = b[1] = c[1] = iy;
          a[2] = b[2] = c[2] = ix;
          for (d = 0; d < 3; d++)
          {
              a[d] += vertList[triTable[cube_idx][t*3]][d];
              b[d] += vertList[triTable[cube_idx][t*3 + 1]][d];
              c[d] += vertList[triTable[cube_idx][t*3 + 2]][d];
              // Factor in the spacing
              a[d] *= spacing[d];
              b[d] *= spacing[d];
              c[d] *= spacing[d];
          }

          // ************************
          // Calculate volume
          // ************************

          // Calculate the cross product
          ab[0] = (a[1] * b[2]) - (b[1] * a[2]);
          ab[1] = (a[2] * b[0]) - (b[2] * a[0]);
          ab[2] = (a[0] * b[1]) - (b[0] * a[1]);

          // Calculate the dot-product and add it to the volume total. The division by 6 is performed at the end.
          *volume += sign_correction * (ab[0] * c[0] + ab[1] * c[1] + ab[2] * c[2]);

          // ************************
          // Calculate surface area
          // ************************

          // Compute the surface, which is equal to 1/2 magnitude of the cross product, where
          // The magnitude is obtained by calculating the euclidean distance between (0, 0, 0)
          // and the location of c
          for (d = 0; d < 3; d++)
          {
            a[d] -= c[d];
            b[d] -= c[d];
          }

          // Compute the cross-product
          ab[0] = (a[1] * b[2]) - (b[1] * a[2]);
          ab[1] = (a[2] * b[0]) - (b[2] * a[0]);
          ab[2] = (a[0] * b[1]) - (b[0] * a[1]);

          // Get the euclidean distance by computing the square and then the square root of the sum.
          ab[0] = ab[0] * ab[0];
          ab[1] = ab[1] * ab[1];
          ab[2] = ab[2] * ab[2];

          sum = ab[0] + ab[1] + ab[2];
          sum = sqrt(sum);

          // multiply by 0.5 (1/2 the magnitude of the cross product)
          sum = 0.5 * sum;

          // Add the surface area of the face to the grand total.
          *surfaceArea += sum;
          t++;
        }
      }
    }
  }
  *volume = *volume / 6;

  // ************************
  // Calculate Diameters using found vertices
  // ************************
  calculate_meshDiameter(vertices, v_idx, diameters);
  free(vertices);
  return 0;
}

void calculate_meshDiameter(double *points, size_t stack_top, double *diameters)
{
  double a[3], b[3], ab[3];
  double distance;
  size_t idx;

  diameters[0] = 0;
  diameters[1] = 0;
  diameters[2] = 0;
  diameters[3] = 0;

  // when the first item is popped, it is the last item entered
  while(stack_top > 0)
  {
    a[2] = points[--stack_top];
    a[1] = points[--stack_top];
    a[0] = points[--stack_top];

    for (idx = 0; idx < stack_top; idx += 3)
    {
      b[0] = points[idx];
      b[1] = points[idx + 1];
      b[2] = points[idx + 2];

      ab[0] = a[0] - b[0];
      ab[1] = a[1] - b[1];
      ab[2] = a[2] - b[2];

      ab[0] *= ab[0];
      ab[1] *= ab[1];
      ab[2] *= ab[2];

      distance = ab[0] + ab[1] + ab[2];
      if (a[0] == b[0] && distance > diameters[0]) diameters[0] = distance;
      if (a[1] == b[1] && distance > diameters[1]) diameters[1] = distance;
      if (a[2] == b[2] && distance > diameters[2]) diameters[2] = distance;
      if (distance > diameters[3]) diameters[3] = distance;
    }
  }

  diameters[0] = sqrt(diameters[0]);
  diameters[1] = sqrt(diameters[1]);
  diameters[2] = sqrt(diameters[2]);
  diameters[3] = sqrt(diameters[3]);
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
  { 1, 9, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 8, 3, 1, 9, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 0, 0, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 1, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 1, 0, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 0, 3, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 8, 10, 11, 10, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 4, 3, 4, 7, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 4, 1, 4, 7, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 10, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 2, 10, 9, 0, 2, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 10, 9, 2, 9, 7, 2, 7, 3, 4, 7, 9, -1, -1, -1, -1 },
  { 4, 7, 8, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 4, 7, 8, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 7, 11, 9, 4, 11, 11, 2, 9, 1, 9, 2, -1, -1, -1, -1 },
  { 3, 10, 1, 11, 10, 3, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 10, 1, 1, 4, 11, 1, 0, 4, 4, 7, 11, -1, -1, -1, -1 },
  { 4, 7, 8, 9, 0, 3, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1 },
  { 4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 4, 0, 1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 3, 1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 10, 5, 5, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 10, 5, 3, 2, 5, 3, 5, 4, 8, 3, 4, -1, -1, -1, -1 },
  { 9, 5, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 0, 0, 8, 11, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 5, 4, 0, 1, 5, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 1, 5, 2, 5, 8, 11, 2, 8, 5, 4, 8, -1, -1, -1, -1 },
  { 3, 10, 1, 11, 10, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 4, 9, 10, 1, 0, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1 },
  { 5, 4, 0, 5, 0, 11, 10, 5, 11, 11, 0, 3, -1, -1, -1, -1 },
  { 5, 4, 8, 10, 5, 8, 11, 10, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 8, 9, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 3, 0, 9, 5, 3, 9, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 8, 0, 0, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 5, 7, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 8, 9, 5, 7, 9, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 10, 3, 0, 9, 3, 9, 5, 3, 5, 7, -1, -1, -1, -1 },
  { 0, 2, 8, 8, 2, 5, 8, 5, 7, 2, 10, 5, -1, -1, -1, -1 },
  { 2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 7, 9, 7, 8, 9, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 7, 9, 9, 7, 2, 2, 0, 9, 11, 2, 7, -1, -1, -1, -1 },
  { 11, 2, 3, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1 },
  { 2, 1, 11, 1, 7, 11, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
  { 7, 8, 9, 5, 7, 9, 3, 10, 1, 11, 10, 3, -1, -1, -1, -1 },
  { 5, 7, 0, 5, 0, 9, 7, 11, 0, 10, 1, 0, 11, 10, 0, -1 },
  { 11, 10, 0, 0, 3, 11, 10, 5, 0, 0, 7, 8, 5, 7, 0, -1 },
  { 11, 10, 5, 5, 7, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 5, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 8, 3, 1, 9, 8, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 6, 5, 1, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 6, 5, 1, 2, 6, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 5, 9, 8, 5, 8, 2, 5, 2, 6, 8, 3, 2, -1, -1, -1, -1 },
  { 11, 2, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 0, 0, 8, 11, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 11, 2, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 1, 1, 9, 11, 9, 8, 11, 6, 5, 10, -1, -1, -1, -1 },
  { 6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
  { 6, 3, 11, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1 },
  { 6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
  { 4, 7, 8, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 4, 3, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 0, 4, 7, 8, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 9, 4, 1, 4, 7, 1, 7, 3, 6, 5, 10, -1, -1, -1, -1 },
  { 4, 7, 8, 1, 6, 5, 1, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 4, 3, 4, 7, 3, 1, 6, 5, 1, 2, 6, -1, -1, -1, -1 },
  { 4, 7, 8, 9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1 },
  { 7, 3, 9, 4, 7, 9, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1 },
  { 11, 2, 3, 4, 7, 8, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 4, 7, 11, 2, 4, 2, 0, 4, 6, 5, 10, -1, -1, -1, -1 },
  { 1, 9, 0, 11, 2, 3, 4, 7, 8, 6, 5, 10, -1, -1, -1, -1 },
  { 4, 7, 11, 9, 4, 11, 11, 2, 9, 1, 9, 2, 6, 5, 10, -1 },
  { 4, 7, 8, 6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1 },
  { 5, 1, 11, 5, 11, 6, 1, 0, 11, 4, 7, 11, 0, 4, 11, -1 },
  { 4, 7, 8, 6, 3, 11, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1 },
  { 6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1 },
  { 10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 10, 4, 9, 6, 4, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 1, 10, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 8, 3, 1, 6, 8, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
  { 1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 4, 9, 1, 2, 4, 2, 6, 4, 0, 8, 3, -1, -1, -1, -1 },
  { 0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 3, 10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 0, 0, 8, 11, 10, 4, 9, 6, 4, 10, -1, -1, -1, -1 },
  { 11, 2, 3, 0, 1, 10, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1 },
  { 6, 4, 1, 6, 1, 10, 1, 4, 8, 11, 2, 1, 1, 8, 11, -1 },
  { 9, 6, 4, 3, 6, 9, 1, 3, 9, 11, 6, 3, -1, -1, -1, -1 },
  { 1, 8, 11, 0, 8, 1, 11, 6, 1, 1, 4, 9, 6, 4, 1, -1 },
  { 6, 3, 11, 0, 3, 6, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
  { 8, 6, 4, 6, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 7, 10, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
  { 6, 7, 10, 0, 10, 7, 0, 9, 10, 0, 7, 3, -1, -1, -1, -1 },
  { 6, 7, 10, 1, 10, 7, 1, 7, 8, 0, 1, 8, -1, -1, -1, -1 },
  { 6, 7, 10, 1, 10, 7, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
  { 1, 2, 6, 1, 6, 8, 1, 8, 9, 6, 7, 8, -1, -1, -1, -1 },
  { 2, 6, 9, 1, 2, 9, 6, 7, 9, 3, 0, 9, 7, 3, 9, -1 },
  { 0, 7, 8, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
  { 2, 7, 3, 2, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 11, 2, 3, 6, 7, 10, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1 },
  { 2, 0, 7, 11, 2, 7, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1 },
  { 6, 7, 10, 1, 10, 7, 1, 7, 8, 0, 1, 8, 11, 2, 3, -1 },
  { 11, 2, 1, 1, 7, 11, 10, 6, 1, 1, 6, 7, -1, -1, -1, -1 },
  { 8, 9, 6, 6, 7, 8, 1, 6, 9, 11, 6, 3, 1, 3, 6, -1 },
  { 0, 9, 1, 6, 7, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
  { 0, 7, 8, 7, 0, 6, 0, 3, 11, 11, 6, 0, -1, -1, -1, -1 },
  { 6, 7, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 }
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


// **************************************************

// 2D Shape calculation

// **************************************************


double calculate_meshDiameter2D(double *points, size_t stack_top);

// Declare the look-up tables, these are filled at the bottom of this code file.
static const int gridAngles2D[4][2];
static const int lineTable2D[16][5];
static const double vertList2D[4][2];

int calculate_coefficients2D(char *mask, int *size, int *strides, double *spacing,
                             double *perimeter, double *surface, double *diameter)
{
  int iy, ix, i, t, d;  // iterator indices
  unsigned char square_idx;  // cube identifier, 4 bits signifying which corners of the cube belong to the segmentation
  int a_idx;  // Angle index (4 'angles', one pointing to each corner of the marching cube

  static const int points_edges[2][2] = {{0, 2}, {3, 2}};
  size_t v_idx = 0;
  size_t v_max = 0;
  double *vertices;

  double sum;
  double a[2], b[2];  // 2 points of the line

  *perimeter = 0;  // Total perimeter
  *surface = 0;  // Total surface

  // create a stack to hold the found vertices. For each square, a maximum of 2 vertices are stored (with x and y)
  // coordinates). This prevents double storing of the vertices.
  v_max = (size[0] - 1) * (size[1] - 1) * 4;
  vertices = (double *)calloc(v_max, sizeof(double));

  // Iterate over all pixels, do not include last voxels in the three dimensions, as the cube includes voxels at pos +1
  for (iy = 0; iy < (size[0] - 1); iy++)
  {
    for (ix = 0; ix < (size[1] - 1); ix++)
    {

      /* Get current square_idx by analyzing each point of the current square (origin is in left-upper corner)
      *  O - X
      *  |
      *  Y
      *         v0
      *   p0 ------- p1
      *    |         |
      * v3 |         | v1
      *    |         |
      *   p3 ------- p2
      *         v2
      */

      square_idx = 0;
      for (a_idx = 0; a_idx < 4; a_idx++)
      {
        i = (iy + gridAngles2D[a_idx][0]) * strides[0] +
            (ix + gridAngles2D[a_idx][1]) * strides[1];

        if (mask[i])
          square_idx |= (1 << a_idx);
      }

      // Exlcude squares entirely outside or inside the segmentation (square_idx = 0 or 0xF = B1111).
      if (square_idx == 0 || square_idx == 0xF)
        continue;

      // Process all lines for this square
      t = 0;
      while (lineTable2D[square_idx][t*2] >= 0) // Exit loop when no more lines are present (element at index = -1)
      {
        a[0] = b[0] = iy;
        a[1] = b[1] = ix;
        for (d = 0; d < 2; d++)
        {
            a[d] += vertList2D[lineTable2D[square_idx][t*2]][d];
            b[d] += vertList2D[lineTable2D[square_idx][t*2 + 1]][d];
            // Factor in the spacing
            a[d] *= spacing[d];
            b[d] *= spacing[d];
        }

        // ************************
        // Calculate Surface
        // ************************

        // Calculate the cross product. Because for both vectors, z = 0, only the last term need be calculated
        // The surface of the triangle is only 1/2 the magnitude of this result, but the division by 2 is done on the
        // final sum.
        *surface += (a[0] * b[1]) - (b[0] * a[1]);

        // ************************
        // Calculate perimeter
        // ************************

        // Compute the euclidean distance between points a and b.
        // Add the result to the grand total, as the perimeter is the sum of
        // all line lengths.
        for (d = 0; d < 2; d++)
        {
          a[d] -= b[d];

          // Get the euclidean distance by computing the square...
          a[d] = a[d] * a[d];
        }

        // ... and then the square root of the sum.
        sum = a[0] + a[1];
        sum = sqrt(sum);

        // Add the length of the line to the grand total.
        *perimeter += sum;
        t++;
      }

      // ************************
      // Store vertices for diameter calculation
      // ************************

      // check if there are vertices on edges 3 and 2
      // Because of the symmetry around the midpoint and the flip if cube_idx > 0xF, the 4th point will never appear
      // as segmented at this point. Therefore, to check if there are vertices on the adjacent edges (3 and 2),
      // one only needs to check if the corresponding points (0 and 2, respectively) are segmented.
      if (v_idx + 9 > v_max) // Overflow!
      {
        free(vertices);
        return 1;
      }

      if (square_idx > 7)
        square_idx = square_idx ^ 0xF;  // Flip the square index

      for (t = 0; t < 2; t++)
      {
        if (square_idx & (1 << points_edges[0][t]))
        {
          vertices[v_idx++] = (((double)iy) + vertList2D[points_edges[1][t]][0]) * spacing[0];
          vertices[v_idx++] = (((double)ix) + vertList2D[points_edges[1][t]][1]) * spacing[1];
        }
      }
    }
  }
  // The surface area of the triangle is 1/2 the magnitude of the cross product.
  *surface = *surface / 2;

  // ************************
  // Calculate Diameters using found vertices
  // ************************
  *diameter = calculate_meshDiameter2D(vertices, v_idx);
  free(vertices);
  return 0;
}


double calculate_meshDiameter2D(double *points, size_t stack_top)
{
  double diameter = 0;
  double a[2], b[2], ab[3];
  double distance;
  size_t idx;

  // so when the first item is popped, it is the last item entered
  while(stack_top > 0)
  {
    a[1] = points[--stack_top];
    a[0] = points[--stack_top];

    for (idx = 0; idx < stack_top; idx += 2)
    {
      b[0] = points[idx];
      b[1] = points[idx + 1];

      ab[0] = a[0] - b[0];
      ab[1] = a[1] - b[1];

      ab[0] *= ab[0];
      ab[1] *= ab[1];

      distance = ab[0] + ab[1];
      if (distance > diameter)
        diameter = distance;
    }
  }

  return sqrt(diameter);
}


static const int gridAngles2D[4][2] = { { 0, 0 }, { 0, 1 }, { 1, 1 }, { 1, 0 } };


static const int lineTable2D[16][5] = {
  { -1, -1, -1, -1, -1},
  {  3,  0, -1, -1, -1},
  {  0,  1, -1, -1, -1},
  {  3,  1, -1, -1, -1},
  {  1,  2, -1, -1, -1},
  {  1,  2,  3,  0, -1},
  {  0,  2, -1, -1, -1},
  {  3,  2, -1, -1, -1},
  {  2,  3, -1, -1, -1},
  {  2,  0, -1, -1, -1},
  {  0,  1,  2,  3, -1},
  {  2,  1, -1, -1, -1},
  {  1,  3, -1, -1, -1},
  {  1,  0, -1, -1, -1},
  {  0,  3, -1, -1, -1},
  { -1, -1, -1, -1, -1},
};

static const double vertList2D[4][2] = { { 0, 0.5 }, { 0.5, 1 }, { 1, 0.5 }, { 0.5, 0 }};
