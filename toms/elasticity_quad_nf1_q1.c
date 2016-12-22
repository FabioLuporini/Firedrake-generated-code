
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
  static const double FE1_C1_D001[1][2]  = {{-1.0, 1.0}};
  static const double FE0[1][4]  = {{0.25, 0.25, 0.25, 0.25}};
  
  const unsigned int nzc10[2] = {0, 2};
  
  const unsigned int nzc9[2] = {0, 3};
  
  const unsigned int nzc6[2] = {8, 11};
  
  const unsigned int nzc11[2] = {0, 1};
  
  const unsigned int nzc7[2] = {8, 10};
  
  const unsigned int nzc3[2] = {4, 6};
  
  const unsigned int nzc4[2] = {4, 5};
  
  const unsigned int nzc8[2] = {8, 9};
  
  const unsigned int nzc2[2] = {4, 7};
  
  double G[45];
  G[0] = 0.5*K[2]*K[7]*W1*det;
  G[1] = 0.5*K[1]*K[2]*W1*det;
  G[2] = 0.5*K[2]*K[4]*W1*det;
  G[3] = 0.5*K[3]*K[8]*W1*det;
  G[4] = 0.5*K[6]*K[8]*W1*det;
  G[5] = 0.5*K[3]*K[5]*W1*det;
  G[6] = 0.5*K[5]*K[6]*W1*det;
  G[7] = 0.5*K[0]*K[5]*W1*det;
  G[8] = 0.5*K[5]*K[7]*W1*det;
  G[9] = 0.5*K[1]*K[5]*W1*det;
  G[10] = 0.5*K[4]*K[5]*W1*det;
  G[11] = W1*det*(0.5*(K[4]*K[7] + K[5]*K[8]) + K[3]*K[6]);
  G[12] = W1*det*(0.5*(K[1]*K[7] + K[2]*K[8]) + K[0]*K[6]);
  G[13] = W1*det*(0.5*(K[1]*K[1] + K[2]*K[2]) + K[0]*K[0]);
  G[14] = 0.5*K[4]*K[6]*W1*det;
  G[15] = W1*det*(0.5*(K[4]*K[4] + K[5]*K[5]) + K[3]*K[3]);
  G[16] = W1*det*(0.5*(K[1]*K[4] + K[2]*K[5]) + K[0]*K[3]);
  G[17] = 0.5*K[6]*K[7]*W1*det;
  G[18] = W1*det*(0.5*(K[6]*K[6] + K[8]*K[8]) + K[7]*K[7]);
  G[19] = 0.5*K[3]*K[7]*W1*det;
  G[20] = 0.5*K[0]*K[7]*W1*det;
  G[21] = 0.5*K[0]*K[8]*W1*det;
  G[22] = 0.5*K[3]*K[4]*W1*det;
  G[23] = 0.5*K[0]*K[4]*W1*det;
  G[24] = 0.5*K[1]*K[6]*W1*det;
  G[25] = 0.5*K[1]*K[3]*W1*det;
  G[26] = 0.5*K[0]*K[1]*W1*det;
  G[27] = W1*det*(0.5*(K[0]*K[6] + K[2]*K[8]) + K[1]*K[7]);
  G[28] = 0.5*K[4]*K[8]*W1*det;
  G[29] = W1*det*(0.5*(K[3]*K[6] + K[4]*K[7]) + K[5]*K[8]);
  G[30] = W1*det*(0.5*(K[6]*K[6] + K[7]*K[7]) + K[8]*K[8]);
  G[31] = W1*det*(0.5*(K[0]*K[6] + K[1]*K[7]) + K[2]*K[8]);
  G[32] = W1*det*(0.5*(K[0]*K[3] + K[1]*K[4]) + K[2]*K[5]);
  G[33] = W1*det*(0.5*(K[3]*K[6] + K[5]*K[8]) + K[4]*K[7]);
  G[34] = 0.5*K[7]*K[8]*W1*det;
  G[35] = 0.5*K[1]*K[8]*W1*det;
  G[36] = 0.5*K[2]*K[6]*W1*det;
  G[37] = 0.5*K[2]*K[3]*W1*det;
  G[38] = 0.5*K[0]*K[2]*W1*det;
  G[39] = W1*det*(0.5*(K[0]*K[3] + K[2]*K[5]) + K[1]*K[4]);
  G[40] = W1*det*(0.5*(K[0]*K[0] + K[2]*K[2]) + K[1]*K[1]);
  G[41] = W1*det*(0.5*(K[7]*K[7] + K[8]*K[8]) + K[6]*K[6]);
  G[42] = W1*det*(0.5*(K[3]*K[3] + K[5]*K[5]) + K[4]*K[4]);
  G[43] = W1*det*(0.5*(K[0]*K[0] + K[1]*K[1]) + K[2]*K[2]);
  G[44] = W1*det*(0.5*(K[3]*K[3] + K[4]*K[4]) + K[5]*K[5]);
  
  double F0  = 0.0;
  
  
  for (int r  = 0; r < 4; r += 1)
  {
    F0 += (w0[r][0] * FE0[0][r]);
    
  }
  
  double Z[45];
  Z[0] = F0*G[0];
  Z[1] = F0*G[1];
  Z[2] = F0*G[2];
  Z[3] = F0*G[3];
  Z[4] = F0*G[4];
  Z[5] = F0*G[5];
  Z[6] = F0*G[6];
  Z[7] = F0*G[7];
  Z[8] = F0*G[8];
  Z[9] = F0*G[9];
  Z[10] = F0*G[10];
  Z[11] = F0*G[11];
  Z[12] = F0*G[12];
  Z[13] = F0*G[13];
  Z[14] = F0*G[14];
  Z[15] = F0*G[15];
  Z[16] = F0*G[16];
  Z[17] = F0*G[17];
  Z[18] = F0*G[18];
  Z[19] = F0*G[19];
  Z[20] = F0*G[20];
  Z[21] = F0*G[21];
  Z[22] = F0*G[22];
  Z[23] = F0*G[23];
  Z[24] = F0*G[24];
  Z[25] = F0*G[25];
  Z[26] = F0*G[26];
  Z[27] = F0*G[27];
  Z[28] = F0*G[28];
  Z[29] = F0*G[29];
  Z[30] = F0*G[30];
  Z[31] = F0*G[31];
  Z[32] = F0*G[32];
  Z[33] = F0*G[33];
  Z[34] = F0*G[34];
  Z[35] = F0*G[35];
  Z[36] = F0*G[36];
  Z[37] = F0*G[37];
  Z[38] = F0*G[38];
  Z[39] = F0*G[39];
  Z[40] = F0*G[40];
  Z[41] = F0*G[41];
  Z[42] = F0*G[42];
  Z[43] = F0*G[43];
  Z[44] = F0*G[44];
  
  
  for (int j  = 0; j < 2; j += 1)
  {
    
    for (int k  = 0; k < 2; k += 1)
    {
      A[nzc10[j]][nzc10[k]] += (Z[15] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc10[j]][nzc11[k]] += (Z[16] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc10[j]][nzc2[k]] += (Z[14] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc10[j]][nzc3[k]] += (Z[22] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc10[j]][nzc4[k]] += (Z[23] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc10[j]][nzc6[k]] += (Z[6] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc10[j]][nzc7[k]] += (Z[5] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc10[j]][nzc8[k]] += (Z[7] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc10[j]][nzc9[k]] += (Z[11] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc10[k]] += (Z[16] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc11[k]] += (Z[13] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc2[k]] += (Z[24] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc3[k]] += (Z[25] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc4[k]] += (Z[26] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc6[k]] += (Z[36] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc7[k]] += (Z[37] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc8[k]] += (Z[38] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc11[j]][nzc9[k]] += (Z[12] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc10[k]] += (Z[14] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc11[k]] += (Z[24] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc2[k]] += (Z[18] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc3[k]] += (Z[33] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc4[k]] += (Z[27] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc6[k]] += (Z[34] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc7[k]] += (Z[28] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc8[k]] += (Z[35] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc2[j]][nzc9[k]] += (Z[17] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc10[k]] += (Z[22] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc11[k]] += (Z[25] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc2[k]] += (Z[33] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc3[k]] += (Z[42] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc4[k]] += (Z[39] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc6[k]] += (Z[8] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc7[k]] += (Z[10] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc8[k]] += (Z[9] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc3[j]][nzc9[k]] += (Z[19] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc10[k]] += (Z[23] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc11[k]] += (Z[26] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc2[k]] += (Z[27] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc3[k]] += (Z[39] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc4[k]] += (Z[40] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc6[k]] += (Z[0] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc7[k]] += (Z[2] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc8[k]] += (Z[1] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc4[j]][nzc9[k]] += (Z[20] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc10[k]] += (Z[6] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc11[k]] += (Z[36] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc2[k]] += (Z[34] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc3[k]] += (Z[8] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc4[k]] += (Z[0] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc6[k]] += (Z[30] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc7[k]] += (Z[29] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc8[k]] += (Z[31] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc6[j]][nzc9[k]] += (Z[4] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc10[k]] += (Z[5] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc11[k]] += (Z[37] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc2[k]] += (Z[28] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc3[k]] += (Z[10] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc4[k]] += (Z[2] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc6[k]] += (Z[29] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc7[k]] += (Z[44] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc8[k]] += (Z[32] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc7[j]][nzc9[k]] += (Z[3] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc10[k]] += (Z[7] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc11[k]] += (Z[38] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc2[k]] += (Z[35] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc3[k]] += (Z[9] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc4[k]] += (Z[1] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc6[k]] += (Z[31] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc7[k]] += (Z[32] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc8[k]] += (Z[43] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc8[j]][nzc9[k]] += (Z[21] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc10[k]] += (Z[11] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc11[k]] += (Z[12] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc2[k]] += (Z[17] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc3[k]] += (Z[19] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc4[k]] += (Z[20] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc6[k]] += (Z[4] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc7[k]] += (Z[3] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc8[k]] += (Z[21] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      A[nzc9[j]][nzc9[k]] += (Z[41] * FE1_C1_D001[0][k] * FE1_C1_D001[0][j]);
      
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
        
        