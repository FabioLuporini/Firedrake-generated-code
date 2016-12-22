
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[4][4] , double** coordinate_dofs , double** w0 )
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
  
  
  static const double W5[5]  = {-0.133333333333333, 0.075, 0.075, 0.075, 0.075};
  static const double FE0_D001[5][2]  = {{-1.0, 1.0}, 
  {-1.0, 1.0}, 
  {-1.0, 1.0}, 
  {-1.0, 1.0}, 
  {-1.0, 1.0}};
  static const double FE0[5][4]  = {{0.25, 0.25, 0.25, 0.25}, 
  {0.166666666666667, 0.5, 0.166666666666667, 0.166666666666667}, 
  {0.166666666666667, 0.166666666666667, 0.5, 0.166666666666667}, 
  {0.166666666666667, 0.166666666666667, 0.166666666666667, 0.5}, 
  {0.5, 0.166666666666667, 0.166666666666667, 0.166666666666667}};
  
  const unsigned int nzc1[2] = {0, 2};
  
  const unsigned int nzc2[2] = {0, 1};
  
  const unsigned int nzc0[2] = {0, 3};
  
  double G[6];
  G[0] = det*(K[0]*K[3] + K[1]*K[4] + K[2]*K[5]);
  G[1] = det*(K[3]*K[6] + K[4]*K[7] + K[5]*K[8]);
  G[2] = det*(K[3]*K[3] + K[4]*K[4] + K[5]*K[5]);
  G[3] = det*(K[0]*K[0] + K[1]*K[1] + K[2]*K[2]);
  G[4] = det*(K[0]*K[6] + K[1]*K[7] + K[2]*K[8]);
  G[5] = det*(K[6]*K[6] + K[7]*K[7] + K[8]*K[8]);
  
  
  for (int ip  = 0; ip < 5; ip += 1)
  {
    double F0  = 0.0;
    
    for (int r  = 0; r < 4; r += 1)
    {
      F0 += (w0[r][0] * FE0[ip][r]);
      
    }
    double Z[7];
    Z[0] = F0*G[0]*W5[ip];
    Z[1] = F0*G[1]*W5[ip];
    Z[2] = F0*G[2]*W5[ip];
    Z[3] = F0*G[3]*W5[ip];
    Z[4] = F0*G[4]*W5[ip];
    Z[5] = F0*G[5]*W5[ip];
    Z[6] = F0*W5[ip]*det;
    
    for (int j  = 0; j < 2; j += 1)
    {
      
      for (int k  = 0; k < 2; k += 1)
      {
        A[nzc0[j]][nzc0[k]] += (Z[5] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        A[nzc0[j]][nzc1[k]] += (Z[1] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        A[nzc0[j]][nzc2[k]] += (Z[4] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        A[nzc1[j]][nzc0[k]] += (Z[1] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        A[nzc1[j]][nzc1[k]] += (Z[2] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        A[nzc1[j]][nzc2[k]] += (Z[0] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        A[nzc2[j]][nzc0[k]] += (Z[4] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        A[nzc2[j]][nzc1[k]] += (Z[0] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        A[nzc2[j]][nzc2[k]] += (Z[3] * FE0_D001[ip][k] * FE0_D001[ip][j]);
        
      }
      
    }
    
    for (int j  = 0; j < 4; j += 1)
    {
      
      for (int k  = 0; k < 4; k += 1)
      {
        A[j][k] += (Z[6] * FE0[ip][k] * FE0[ip][j]);
        
      }
      
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
    double buffer_arg0_0[4][4] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec, arg2_0_vec);
    MatSetValuesLocal(arg0_0_0, 4, arg0_0_map0_0 + i * 4,
                                             4, arg0_0_map1_0 + i * 4,
                                             (const PetscScalar *)buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        