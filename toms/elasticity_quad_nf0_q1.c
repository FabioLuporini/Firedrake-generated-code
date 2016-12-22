
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[12][12] , double** coordinate_dofs )
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
  static const double FE0_C1_D001[1][2]  = {{-1.0, 1.0}};
  
  const unsigned int nzc6[2] = {8, 11};
  
  const unsigned int nzc2[2] = {4, 7};
  
  const unsigned int nzc3[2] = {4, 6};
  
  const unsigned int nzc7[2] = {8, 10};
  
  const unsigned int nzc4[2] = {4, 5};
  
  const unsigned int nzc11[2] = {0, 1};
  
  const unsigned int nzc9[2] = {0, 3};
  
  const unsigned int nzc10[2] = {0, 2};
  
  const unsigned int nzc8[2] = {8, 9};
  
  double G[45];
  G[0] = W1*det*(0.5*(K[1]*K[1] + K[2]*K[2]) + K[0]*K[0]);
  G[1] = W1*det*(0.5*(K[1]*K[4] + K[2]*K[5]) + K[0]*K[3]);
  G[2] = W1*det*(0.5*(K[1]*K[7] + K[2]*K[8]) + K[0]*K[6]);
  G[3] = W1*det*(0.5*(K[4]*K[4] + K[5]*K[5]) + K[3]*K[3]);
  G[4] = W1*det*(0.5*(K[4]*K[7] + K[5]*K[8]) + K[3]*K[6]);
  G[5] = W1*det*(0.5*(K[7]*K[7] + K[8]*K[8]) + K[6]*K[6]);
  G[6] = 0.5*K[0]*K[1]*W1*det;
  G[7] = 0.5*K[1]*K[3]*W1*det;
  G[8] = 0.5*K[1]*K[6]*W1*det;
  G[9] = 0.5*K[0]*K[4]*W1*det;
  G[10] = 0.5*K[3]*K[4]*W1*det;
  G[11] = 0.5*K[4]*K[6]*W1*det;
  G[12] = 0.5*K[0]*K[7]*W1*det;
  G[13] = 0.5*K[3]*K[7]*W1*det;
  G[14] = 0.5*K[6]*K[7]*W1*det;
  G[15] = W1*det*(0.5*(K[0]*K[0] + K[2]*K[2]) + K[1]*K[1]);
  G[16] = W1*det*(0.5*(K[0]*K[3] + K[2]*K[5]) + K[1]*K[4]);
  G[17] = W1*det*(0.5*(K[0]*K[6] + K[2]*K[8]) + K[1]*K[7]);
  G[18] = W1*det*(0.5*(K[3]*K[3] + K[5]*K[5]) + K[4]*K[4]);
  G[19] = W1*det*(0.5*(K[3]*K[6] + K[5]*K[8]) + K[4]*K[7]);
  G[20] = W1*det*(0.5*(K[6]*K[6] + K[8]*K[8]) + K[7]*K[7]);
  G[21] = 0.5*K[0]*K[2]*W1*det;
  G[22] = 0.5*K[2]*K[3]*W1*det;
  G[23] = 0.5*K[2]*K[6]*W1*det;
  G[24] = 0.5*K[0]*K[5]*W1*det;
  G[25] = 0.5*K[3]*K[5]*W1*det;
  G[26] = 0.5*K[5]*K[6]*W1*det;
  G[27] = 0.5*K[0]*K[8]*W1*det;
  G[28] = 0.5*K[3]*K[8]*W1*det;
  G[29] = 0.5*K[6]*K[8]*W1*det;
  G[30] = W1*det*(0.5*(K[0]*K[0] + K[1]*K[1]) + K[2]*K[2]);
  G[31] = W1*det*(0.5*(K[0]*K[3] + K[1]*K[4]) + K[2]*K[5]);
  G[32] = W1*det*(0.5*(K[0]*K[6] + K[1]*K[7]) + K[2]*K[8]);
  G[33] = W1*det*(0.5*(K[3]*K[3] + K[4]*K[4]) + K[5]*K[5]);
  G[34] = W1*det*(0.5*(K[3]*K[6] + K[4]*K[7]) + K[5]*K[8]);
  G[35] = W1*det*(0.5*(K[6]*K[6] + K[7]*K[7]) + K[8]*K[8]);
  G[36] = 0.5*K[1]*K[8]*W1*det;
  G[37] = 0.5*K[1]*K[5]*W1*det;
  G[38] = 0.5*K[1]*K[2]*W1*det;
  G[39] = 0.5*K[2]*K[4]*W1*det;
  G[40] = 0.5*K[2]*K[7]*W1*det;
  G[41] = 0.5*K[4]*K[8]*W1*det;
  G[42] = 0.5*K[4]*K[5]*W1*det;
  G[43] = 0.5*K[7]*K[8]*W1*det;
  G[44] = 0.5*K[5]*K[7]*W1*det;
  
  
  for (int j  = 0; j < 2; j += 1)
  {
    
    for (int k  = 0; k < 2; k += 1)
    {
      A[nzc10[j]][nzc10[k]] += (G[3] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc10[j]][nzc11[k]] += (G[1] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc10[j]][nzc2[k]] += (G[11] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc10[j]][nzc3[k]] += (G[10] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc10[j]][nzc4[k]] += (G[9] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc10[j]][nzc6[k]] += (G[26] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc10[j]][nzc7[k]] += (G[25] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc10[j]][nzc8[k]] += (G[24] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc10[j]][nzc9[k]] += (G[4] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc10[k]] += (G[1] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc11[k]] += (G[0] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc2[k]] += (G[8] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc3[k]] += (G[7] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc4[k]] += (G[6] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc6[k]] += (G[23] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc7[k]] += (G[22] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc8[k]] += (G[21] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc11[j]][nzc9[k]] += (G[2] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc10[k]] += (G[11] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc11[k]] += (G[8] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc2[k]] += (G[20] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc3[k]] += (G[19] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc4[k]] += (G[17] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc6[k]] += (G[43] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc7[k]] += (G[41] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc8[k]] += (G[36] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc2[j]][nzc9[k]] += (G[14] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc10[k]] += (G[10] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc11[k]] += (G[7] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc2[k]] += (G[19] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc3[k]] += (G[18] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc4[k]] += (G[16] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc6[k]] += (G[44] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc7[k]] += (G[42] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc8[k]] += (G[37] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc3[j]][nzc9[k]] += (G[13] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc10[k]] += (G[9] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc11[k]] += (G[6] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc2[k]] += (G[17] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc3[k]] += (G[16] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc4[k]] += (G[15] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc6[k]] += (G[40] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc7[k]] += (G[39] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc8[k]] += (G[38] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc4[j]][nzc9[k]] += (G[12] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc10[k]] += (G[26] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc11[k]] += (G[23] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc2[k]] += (G[43] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc3[k]] += (G[44] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc4[k]] += (G[40] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc6[k]] += (G[35] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc7[k]] += (G[34] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc8[k]] += (G[32] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc6[j]][nzc9[k]] += (G[29] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc10[k]] += (G[25] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc11[k]] += (G[22] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc2[k]] += (G[41] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc3[k]] += (G[42] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc4[k]] += (G[39] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc6[k]] += (G[34] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc7[k]] += (G[33] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc8[k]] += (G[31] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc7[j]][nzc9[k]] += (G[28] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc10[k]] += (G[24] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc11[k]] += (G[21] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc2[k]] += (G[36] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc3[k]] += (G[37] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc4[k]] += (G[38] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc6[k]] += (G[32] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc7[k]] += (G[31] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc8[k]] += (G[30] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc8[j]][nzc9[k]] += (G[27] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc10[k]] += (G[4] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc11[k]] += (G[2] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc2[k]] += (G[14] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc3[k]] += (G[13] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc4[k]] += (G[12] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc6[k]] += (G[29] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc7[k]] += (G[28] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc8[k]] += (G[27] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      A[nzc9[j]][nzc9[k]] += (G[5] * FE0_C1_D001[0][k] * FE0_C1_D001[0][j]);
      
    }
    
  }
  
}
            

        
        void wrap_form_cell_integral_0_otherwise(int start, int end,
                      Mat arg0_0_, int *arg0_0_map0_0, int *arg0_0_map1_0, double *arg1_0, int *arg1_0_map0_0
                      ) {
  Mat arg0_0_0 = arg0_0_;
  double *arg1_0_vec[12];
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
    double buffer_arg0_0[12][12] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec);
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
        
        