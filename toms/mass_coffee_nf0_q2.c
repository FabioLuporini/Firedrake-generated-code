
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[10][10] , double** coordinate_dofs )
{
  // This code is generated visiting a COFFEE AST
  
  // Preevaluated tables
  static const double IP_J_K_1_2_0[10][10]  = {{0.00238095238095, 0.000396825396825, 0.000396825396825, 0.000396825396825, -0.00238095238095, -0.00238095238095, -0.00238095238095, -0.0015873015873, -0.0015873015873, -0.0015873015873}, 
  {0.000396825396825, 0.00238095238095, 0.000396825396825, 0.000396825396825, -0.00238095238095, -0.0015873015873, -0.0015873015873, -0.00238095238095, -0.00238095238095, -0.0015873015873}, 
  {0.000396825396825, 0.000396825396825, 0.00238095238095, 0.000396825396825, -0.0015873015873, -0.00238095238095, -0.0015873015873, -0.00238095238095, -0.0015873015873, -0.00238095238095}, 
  {0.000396825396825, 0.000396825396825, 0.000396825396825, 0.00238095238095, -0.0015873015873, -0.0015873015873, -0.00238095238095, -0.0015873015873, -0.00238095238095, -0.00238095238095}, 
  {-0.00238095238095, -0.00238095238095, -0.0015873015873, -0.0015873015873, 0.0126984126984, 0.00634920634921, 0.00634920634921, 0.00634920634921, 0.00634920634921, 0.0031746031746}, 
  {-0.00238095238095, -0.0015873015873, -0.00238095238095, -0.0015873015873, 0.00634920634921, 0.0126984126984, 0.00634920634921, 0.00634920634921, 0.0031746031746, 0.00634920634921}, 
  {-0.00238095238095, -0.0015873015873, -0.0015873015873, -0.00238095238095, 0.00634920634921, 0.00634920634921, 0.0126984126984, 0.0031746031746, 0.00634920634921, 0.00634920634921}, 
  {-0.0015873015873, -0.00238095238095, -0.00238095238095, -0.0015873015873, 0.00634920634921, 0.00634920634921, 0.0031746031746, 0.0126984126984, 0.00634920634921, 0.00634920634921}, 
  {-0.0015873015873, -0.00238095238095, -0.0015873015873, -0.00238095238095, 0.00634920634921, 0.0031746031746, 0.00634920634921, 0.00634920634921, 0.0126984126984, 0.00634920634921}, 
  {-0.0015873015873, -0.0015873015873, -0.00238095238095, -0.00238095238095, 0.0031746031746, 0.00634920634921, 0.00634920634921, 0.00634920634921, 0.00634920634921, 0.0126984126984}};
  
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
  
  
  static const double W14[14]  = {0.00317460317460317, 0.00317460317460317, 0.00317460317460317, 0.00317460317460317, 0.00317460317460317, 0.00317460317460317, 0.0147649707904968, 0.0147649707904968, 0.0147649707904968, 0.0147649707904968, 0.0221397911142651, 0.0221397911142651, 0.0221397911142651, 0.0221397911142651};
  static const double FE0[14][10]  = {{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0}, 
  {-0.0803155041719177, 0.277160462452741, -0.0803155041719177, -0.0803155041719177, 0.0404225221065735, 0.280839494581097, 0.280839494581097, 0.0404225221065735, 0.0404225221065734, 0.280839494581097}, 
  {0.27716046245274, -0.0803155041719176, -0.0803155041719177, -0.0803155041719177, 0.0404225221065735, 0.0404225221065736, 0.0404225221065735, 0.280839494581097, 0.280839494581097, 0.280839494581097}, 
  {-0.0803155041719177, -0.0803155041719177, -0.0803155041719177, 0.277160462452741, 0.280839494581097, 0.280839494581097, 0.0404225221065736, 0.280839494581097, 0.0404225221065734, 0.0404225221065735}, 
  {-0.0803155041719177, -0.0803155041719177, 0.277160462452741, -0.0803155041719177, 0.280839494581097, 0.0404225221065736, 0.280839494581097, 0.0404225221065735, 0.280839494581097, 0.0404225221065735}, 
  {-0.116712266316459, -0.0504103968481305, -0.116712266316459, -0.116712266316459, 0.395321214353467, 0.0715278509123693, 0.0715278509123693, 0.395321214353467, 0.395321214353466, 0.0715278509123692}, 
  {-0.0504103968481305, -0.116712266316459, -0.116712266316459, -0.116712266316459, 0.395321214353467, 0.395321214353467, 0.395321214353467, 0.0715278509123693, 0.0715278509123691, 0.0715278509123692}, 
  {-0.116712266316459, -0.116712266316459, -0.116712266316459, -0.0504103968481305, 0.0715278509123693, 0.0715278509123693, 0.395321214353466, 0.0715278509123693, 0.395321214353466, 0.395321214353466}, 
  {-0.116712266316459, -0.116712266316459, -0.0504103968481305, -0.116712266316459, 0.0715278509123693, 0.395321214353467, 0.0715278509123693, 0.395321214353466, 0.0715278509123692, 0.395321214353466}};
  
  
  for (int j  = 0; j < 10; j += 1)
  {
    
    for (int k  = 0; k < 10; k += 1)
    {
      A[j][k] += (det * IP_J_K_1_2_0[j][k]);
      
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
    double buffer_arg0_0[10][10]  = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec);
    MatSetValuesLocal(arg0_0_0, 10, arg0_0_map0_0 + i * 10,
                                             10, arg0_0_map1_0 + i * 10,
                                             (const PetscScalar *)buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        