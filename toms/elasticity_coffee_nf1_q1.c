
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[12][12] , double** coordinate_dofs , double** w0 )
{
  // This code is generated visiting a COFFEE AST
  
  // Compute Jacobian
  double J[9];
  compute_jacobian_tetrahedron_3d(J, coordinate_dofs);
  
  // Compute Jacobian inverse and determinant
  double K[9];
  double detJ;
  compute_jacobian_inverse_tetrahedron_3d(K, detJ, J);
  
  const double det = fabs(detJ);
  
  // Compute cell volume
  
  
  // Compute circumradius
  
  
  static const double W1  = 0.166666666666667;
  static const double FE1_D010[1][12]  = {{-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE1_D001[1][12]  = {{-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE1_C2_D001[1][12]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0}};
  static const double FE1_D100[1][12]  = {{-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE1_C2_D010[1][12]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0}};
  static const double FE1_C1_D010[1][12]  = {{0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE1_C1_D100[1][12]  = {{0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE1_C2_D100[1][12]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0}};
  static const double FE0[1][4]  = {{0.25, 0.25, 0.25, 0.25}};
  static const double FE1_C1_D001[1][12]  = {{0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}};
  double F0  = 0.0;
  
  
  for (int r  = 0; r < 4; r += 1)
  {
    F0 += (w0[r][0] * FE0[0][r]);
    
  }
  
  double CONST_0_1_0  = (det * W1 * F0 * 0.25);
  double J_0_3_0[12]  = {0.0};
  double J_0_3_1[12]  = {0.0};
  double J_0_3_2[12]  = {0.0};
  double J_0_3_3[12]  = {0.0};
  double J_0_3_4[12]  = {0.0};
  double J_0_3_5[12]  = {0.0};
  double K_0_3_0[12]  = {0.0};
  double K_0_3_1[12]  = {0.0};
  double K_0_3_2[12]  = {0.0};
  double K_0_3_3[12]  = {0.0};
  double K_0_3_4[12]  = {0.0};
  double K_0_3_5[12]  = {0.0};
  
  
  for (int k  = 0; k < 4; k += 1)
  {
    J_0_3_0[k+4] = (K[1] * FE1_C1_D100[0][k+4]) + (K[4] * FE1_C1_D010[0][k+4]) + (K[7] * FE1_C1_D001[0][k+4]);
    J_0_3_1[k] = (K[2] * FE1_D100[0][k]) + (K[5] * FE1_D010[0][k]) + (K[8] * FE1_D001[0][k]) + (K[0] * FE1_C2_D100[0][k]) + (K[3] * FE1_C2_D010[0][k]) + (K[6] * FE1_C2_D001[0][k]);
    J_0_3_1[k+8] = (K[2] * FE1_D100[0][k+8]) + (K[5] * FE1_D010[0][k+8]) + (K[8] * FE1_D001[0][k+8]) + (K[0] * FE1_C2_D100[0][k+8]) + (K[3] * FE1_C2_D010[0][k+8]) + (K[6] * FE1_C2_D001[0][k+8]);
    J_0_3_2[k] = (K[0] * FE1_D100[0][k]) + (K[3] * FE1_D010[0][k]) + (K[6] * FE1_D001[0][k]);
    J_0_3_3[k+8] = (K[2] * FE1_C2_D100[0][k+8]) + (K[5] * FE1_C2_D010[0][k+8]) + (K[8] * FE1_C2_D001[0][k+8]);
    K_0_3_0[k+4] = ((FE1_C1_D100[0][k+4] * 2 * K[1]) + (FE1_C1_D010[0][k+4] * 2 * K[4]) + (FE1_C1_D001[0][k+4] * 2 * K[7])) * 2 * CONST_0_1_0;
    K_0_3_3[k] = ((FE1_D100[0][k] * 2 * K[0]) + (FE1_D010[0][k] * 2 * K[3]) + (FE1_D001[0][k] * 2 * K[6])) * 2 * CONST_0_1_0;
    K_0_3_4[k] = ((FE1_D100[0][k] * K[2]) + (FE1_D010[0][k] * K[5]) + (FE1_D001[0][k] * K[8]) + (FE1_C2_D100[0][k] * K[0]) + (FE1_C2_D010[0][k] * K[3]) + (FE1_C2_D001[0][k] * K[6])) * 2 * CONST_0_1_0;
    K_0_3_4[k+8] = ((FE1_D100[0][k+8] * K[2]) + (FE1_D010[0][k+8] * K[5]) + (FE1_D001[0][k+8] * K[8]) + (FE1_C2_D100[0][k+8] * K[0]) + (FE1_C2_D010[0][k+8] * K[3]) + (FE1_C2_D001[0][k+8] * K[6])) * 2 * CONST_0_1_0;
    K_0_3_5[k+8] = ((FE1_C2_D100[0][k+8] * 2 * K[2]) + (FE1_C2_D010[0][k+8] * 2 * K[5]) + (FE1_C2_D001[0][k+8] * 2 * K[8])) * 2 * CONST_0_1_0;
    
  }
  
  
  for (int k  = 0; k < 8; k += 1)
  {
    J_0_3_4[k+4] = (K[1] * FE1_C2_D100[0][k+4]) + (K[4] * FE1_C2_D010[0][k+4]) + (K[7] * FE1_C2_D001[0][k+4]) + (K[2] * FE1_C1_D100[0][k+4]) + (K[5] * FE1_C1_D010[0][k+4]) + (K[8] * FE1_C1_D001[0][k+4]);
    J_0_3_5[k] = (K[1] * FE1_D100[0][k]) + (K[4] * FE1_D010[0][k]) + (K[7] * FE1_D001[0][k]) + (K[0] * FE1_C1_D100[0][k]) + (K[3] * FE1_C1_D010[0][k]) + (K[6] * FE1_C1_D001[0][k]);
    K_0_3_1[k+4] = ((FE1_C2_D100[0][k+4] * K[1]) + (FE1_C2_D010[0][k+4] * K[4]) + (FE1_C2_D001[0][k+4] * K[7]) + (FE1_C1_D100[0][k+4] * K[2]) + (FE1_C1_D010[0][k+4] * K[5]) + (FE1_C1_D001[0][k+4] * K[8])) * 2 * CONST_0_1_0;
    K_0_3_2[k] = ((FE1_D100[0][k] * K[1]) + (FE1_D010[0][k] * K[4]) + (FE1_D001[0][k] * K[7]) + (FE1_C1_D100[0][k] * K[0]) + (FE1_C1_D010[0][k] * K[3]) + (FE1_C1_D001[0][k] * K[6])) * 2 * CONST_0_1_0;
    
  }
  
  
  
  for (int j  = 0; j < 4; j += 1)
  {
    
    for (int k  = 0; k < 4; k += 1)
    {
      A[j][k] += ((K_0_3_3[k] * J_0_3_2[j])) + ((K_0_3_4[k] * J_0_3_1[j]));
      A[j+4][k+4] += ((K_0_3_0[k+4] * J_0_3_0[j+4]));
      A[j+8][k] += ((K_0_3_4[k] * J_0_3_1[j+8]));
      A[j][k+8] += ((K_0_3_4[k+8] * J_0_3_1[j]));
      A[j+8][k+8] += ((K_0_3_4[k+8] * J_0_3_1[j+8])) + ((K_0_3_5[k+8] * J_0_3_3[j+8]));
      
    }
    
  }
  
  
  for (int j  = 0; j < 8; j += 1)
  {
    
    for (int k  = 0; k < 8; k += 1)
    {
      A[j][k] += ((K_0_3_2[k] * J_0_3_5[j]));
      A[j+4][k+4] += ((K_0_3_1[k+4] * J_0_3_4[j+4]));
      
    }
    
  }
  
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
    double buffer_arg0_0[12][12] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec, arg2_0_vec);
                    double tmp_buffer_arg0_0[12][12] __attribute__((aligned(32)));
                    for ( int j = 0; j < 4; j++ ) {
                       for ( int k = 0; k < 3; k++ ) {
                          for ( int l = 0; l < 4; l++ ) {
                             for ( int m = 0; m < 3; m++ ) {
                                tmp_buffer_arg0_0[3*j + k][3*l + m] = buffer_arg0_0[j + 4*k][l + 4*m];
                             }
                          }
                       }
                    }
    MatSetValuesBlockedLocal(arg0_0_0, 4, arg0_0_map0_0 + i * 4,
                                             4, arg0_0_map1_0 + i * 4,
                                             (const PetscScalar *)tmp_buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        