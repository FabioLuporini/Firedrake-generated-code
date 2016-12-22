
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


void form_cell_integral_0_otherwise(double A[4][4], double **coordinate_dofs, double **w0)
{
  // Number of operations (multiply-add pairs) for Jacobian data:      2
  // Number of operations (multiply-add pairs) for geometry tensor:    4
  // Number of operations (multiply-add pairs) for tensor contraction: 56
  // Total number of operations (multiply-add pairs):                  62
  
  // Compute Jacobian
  double J[9];
  compute_jacobian_tetrahedron_3d(J, coordinate_dofs);
  
  // Compute Jacobian inverse and determinant
  double K[9];
  double detJ;
  compute_jacobian_inverse_tetrahedron_3d(K, detJ, J);
  
  const double det = fabs(detJ);
  
  // Compute geometry tensor
  const double G0_0 = det*w0[0][0]*(1.0);
  const double G0_1 = det*w0[1][0]*(1.0);
  const double G0_2 = det*w0[2][0]*(1.0);
  const double G0_3 = det*w0[3][0]*(1.0);
  
  // Compute element tensor
  A[0][0] = 0.00833333333333334*G0_0 + 0.00277777777777778*G0_1 + 0.00277777777777778*G0_2 + 0.00277777777777778*G0_3;
  A[0][1] = 0.00277777777777778*G0_0 + 0.00277777777777778*G0_1 + 0.00138888888888889*G0_2 + 0.00138888888888889*G0_3;
  A[0][2] = 0.00277777777777778*G0_0 + 0.00138888888888889*G0_1 + 0.00277777777777778*G0_2 + 0.00138888888888889*G0_3;
  A[0][3] = 0.00277777777777778*G0_0 + 0.00138888888888889*G0_1 + 0.00138888888888889*G0_2 + 0.00277777777777778*G0_3;
  A[1][0] = 0.00277777777777778*G0_0 + 0.00277777777777778*G0_1 + 0.00138888888888889*G0_2 + 0.00138888888888889*G0_3;
  A[1][1] = 0.00277777777777778*G0_0 + 0.00833333333333333*G0_1 + 0.00277777777777777*G0_2 + 0.00277777777777777*G0_3;
  A[1][2] = 0.00138888888888889*G0_0 + 0.00277777777777778*G0_1 + 0.00277777777777778*G0_2 + 0.00138888888888889*G0_3;
  A[1][3] = 0.00138888888888889*G0_0 + 0.00277777777777777*G0_1 + 0.00138888888888889*G0_2 + 0.00277777777777778*G0_3;
  A[2][0] = 0.00277777777777778*G0_0 + 0.00138888888888889*G0_1 + 0.00277777777777778*G0_2 + 0.00138888888888889*G0_3;
  A[2][1] = 0.00138888888888889*G0_0 + 0.00277777777777778*G0_1 + 0.00277777777777778*G0_2 + 0.00138888888888889*G0_3;
  A[2][2] = 0.00277777777777778*G0_0 + 0.00277777777777778*G0_1 + 0.00833333333333333*G0_2 + 0.00277777777777778*G0_3;
  A[2][3] = 0.00138888888888889*G0_0 + 0.00138888888888889*G0_1 + 0.00277777777777778*G0_2 + 0.00277777777777778*G0_3;
  A[3][0] = 0.00277777777777778*G0_0 + 0.00138888888888889*G0_1 + 0.00138888888888889*G0_2 + 0.00277777777777778*G0_3;
  A[3][1] = 0.00138888888888889*G0_0 + 0.00277777777777777*G0_1 + 0.00138888888888889*G0_2 + 0.00277777777777778*G0_3;
  A[3][2] = 0.00138888888888889*G0_0 + 0.00138888888888889*G0_1 + 0.00277777777777778*G0_2 + 0.00277777777777778*G0_3;
  A[3][3] = 0.00277777777777778*G0_0 + 0.00277777777777778*G0_1 + 0.00277777777777778*G0_2 + 0.00833333333333333*G0_3;
}

            

        
        void wrap_form_cell_integral_0_otherwise(int start, int end,
                      Mat arg0_0_, int *arg0_0_map0_0, int *arg0_0_map1_0, double *arg1_0, int *arg1_0_map0_0, double *arg2_0, int *arg2_0_map0_0
                      ) {
  Mat arg0_0_0 = arg0_0_;
  double *arg1_0_vec[12];
    double *arg2_0_vec[4];
  for ( int n = start; n < end; n++ ) {
    int i = n;
    arg1_0_vec[0] = arg1_0 + (arg1_0_map0_0[i * 4 + 0])* 3;
    arg1_0_vec[1] = arg1_0 + (arg1_0_map0_0[i * 4 + 1])* 3;
    arg1_0_vec[2] = arg1_0 + (arg1_0_map0_0[i * 4 + 2])* 3;
    arg1_0_vec[3] = arg1_0 + (arg1_0_map0_0[i * 4 + 3])* 3;
    arg1_0_vec[4] = arg1_0 + (arg1_0_map0_0[i * 4 + 0])* 3 + 1;
    arg1_0_vec[5] = arg1_0 + (arg1_0_map0_0[i * 4 + 1])* 3 + 1;
    arg1_0_vec[6] = arg1_0 + (arg1_0_map0_0[i * 4 + 2])* 3 + 1;
    arg1_0_vec[7] = arg1_0 + (arg1_0_map0_0[i * 4 + 3])* 3 + 1;
    arg1_0_vec[8] = arg1_0 + (arg1_0_map0_0[i * 4 + 0])* 3 + 2;
    arg1_0_vec[9] = arg1_0 + (arg1_0_map0_0[i * 4 + 1])* 3 + 2;
    arg1_0_vec[10] = arg1_0 + (arg1_0_map0_0[i * 4 + 2])* 3 + 2;
    arg1_0_vec[11] = arg1_0 + (arg1_0_map0_0[i * 4 + 3])* 3 + 2;
    arg2_0_vec[0] = arg2_0 + (arg2_0_map0_0[i * 4 + 0])* 1;
    arg2_0_vec[1] = arg2_0 + (arg2_0_map0_0[i * 4 + 1])* 1;
    arg2_0_vec[2] = arg2_0 + (arg2_0_map0_0[i * 4 + 2])* 1;
    arg2_0_vec[3] = arg2_0 + (arg2_0_map0_0[i * 4 + 3])* 1;
    double buffer_arg0_0[4][4] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec, arg2_0_vec);
    MatSetValuesLocal(arg0_0_0, 4, arg0_0_map0_0 + i * 4,
                                             4, arg0_0_map1_0 + i * 4,
                                             (const PetscScalar *)buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        