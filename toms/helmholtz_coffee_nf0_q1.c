
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[4][4] , double** coordinate_dofs )
{
  // This code is generated visiting a COFFEE AST
  
  // Preevaluated tables
  static const double IP_J_K_1_3_9[4][4]  = {{0.0166666666667, 0.00833333333333, 0.00833333333333, 0.00833333333333}, 
  {0.00833333333333, 0.0166666666667, 0.00833333333333, 0.00833333333333}, 
  {0.00833333333333, 0.00833333333333, 0.0166666666667, 0.00833333333333}, 
  {0.00833333333333, 0.00833333333333, 0.00833333333333, 0.0166666666667}};
  static const double IP_J_K_1_3_8[4][4]  = {{0.166666666667, 0.0, 0.0, -0.166666666667}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {-0.166666666667, 0.0, 0.0, 0.166666666667}};
  static const double IP_J_K_1_3_7[4][4]  = {{0.166666666667, 0.0, 0.0, -0.166666666667}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {-0.166666666667, 0.0, 0.0, 0.166666666667}, 
  {0.0, 0.0, 0.0, 0.0}};
  static const double IP_J_K_1_3_6[4][4]  = {{0.166666666667, 0.0, -0.166666666667, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {-0.166666666667, 0.0, 0.166666666667, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}};
  static const double IP_J_K_1_3_5[4][4]  = {{0.166666666667, 0.0, 0.0, -0.166666666667}, 
  {-0.166666666667, 0.0, 0.0, 0.166666666667}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}};
  static const double IP_J_K_1_3_4[4][4]  = {{0.166666666667, 0.0, -0.166666666667, 0.0}, 
  {-0.166666666667, 0.0, 0.166666666667, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}};
  static const double IP_J_K_1_3_3[4][4]  = {{0.166666666667, -0.166666666667, 0.0, 0.0}, 
  {-0.166666666667, 0.166666666667, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}};
  static const double IP_J_K_1_3_2[4][4]  = {{0.166666666667, 0.0, -0.166666666667, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {-0.166666666667, 0.0, 0.166666666667, 0.0}};
  static const double IP_J_K_1_3_1[4][4]  = {{0.166666666667, -0.166666666667, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {-0.166666666667, 0.166666666667, 0.0, 0.0}};
  static const double IP_J_K_1_3_0[4][4]  = {{0.166666666667, -0.166666666667, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}, 
  {-0.166666666667, 0.166666666667, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0}};
  
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
  
  
  static const double W4[4]  = {0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 0.0416666666666667};
  static const double FE0_D100[4][4]  = {{-1.0, 1.0, 0.0, 0.0}, 
  {-1.0, 1.0, 0.0, 0.0}, 
  {-1.0, 1.0, 0.0, 0.0}, 
  {-1.0, 1.0, 0.0, 0.0}};
  static const double FE0_D001[4][4]  = {{-1.0, 0.0, 0.0, 1.0}, 
  {-1.0, 0.0, 0.0, 1.0}, 
  {-1.0, 0.0, 0.0, 1.0}, 
  {-1.0, 0.0, 0.0, 1.0}};
  static const double FE0[4][4]  = {{0.138196601125009, 0.585410196624969, 0.138196601125011, 0.138196601125011}, 
  {0.138196601125009, 0.138196601125011, 0.585410196624969, 0.138196601125011}, 
  {0.138196601125009, 0.138196601125011, 0.138196601125011, 0.585410196624969}, 
  {0.585410196624967, 0.138196601125011, 0.138196601125011, 0.138196601125011}};
  static const double FE0_D010[4][4]  = {{-1.0, 0.0, 1.0, 0.0}, 
  {-1.0, 0.0, 1.0, 0.0}, 
  {-1.0, 0.0, 1.0, 0.0}, 
  {-1.0, 0.0, 1.0, 0.0}};
  double CONST_1_1_0  = ((K[8] * K[5] * det) + (K[7] * K[4] * det) + (K[6] * K[3] * det));
  double CONST_1_1_1  = ((K[2] * K[8] * det) + (K[1] * K[7] * det) + (K[0] * K[6] * det));
  double CONST_1_1_2  = ((K[2] * K[2] * det) + (K[1] * K[1] * det) + (K[0] * K[0] * det));
  double CONST_1_1_3  = ((K[5] * K[2] * det) + (K[4] * K[1] * det) + (K[3] * K[0] * det));
  double CONST_1_1_4  = ((K[2] * K[5] * det) + (K[1] * K[4] * det) + (K[0] * K[3] * det));
  double CONST_1_1_5  = ((K[8] * K[2] * det) + (K[7] * K[1] * det) + (K[6] * K[0] * det));
  double CONST_1_1_6  = ((K[5] * K[5] * det) + (K[4] * K[4] * det) + (K[3] * K[3] * det));
  double CONST_1_1_7  = ((K[8] * K[8] * det) + (K[7] * K[7] * det) + (K[6] * K[6] * det));
  double CONST_1_1_8  = ((K[5] * K[8] * det) + (K[4] * K[7] * det) + (K[3] * K[6] * det));
  
  
  for (int j  = 0; j < 4; j += 1)
  {
    
    for (int k  = 0; k < 4; k += 1)
    {
      A[j][k] += (((IP_J_K_1_3_9[j][k] * det) + (IP_J_K_1_3_3[j][k] * CONST_1_1_2) + (IP_J_K_1_3_0[j][k] * CONST_1_1_4) + (IP_J_K_1_3_1[j][k] * CONST_1_1_1) + (IP_J_K_1_3_4[j][k] * CONST_1_1_3) + (IP_J_K_1_3_6[j][k] * CONST_1_1_6) + (IP_J_K_1_3_2[j][k] * CONST_1_1_8) + (IP_J_K_1_3_5[j][k] * CONST_1_1_5) + (IP_J_K_1_3_7[j][k] * CONST_1_1_0) + (IP_J_K_1_3_8[j][k] * CONST_1_1_7)));
      
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
    double buffer_arg0_0[4][4] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec);
    MatSetValuesLocal(arg0_0_0, 4, arg0_0_map0_0 + i * 4,
                                             4, arg0_0_map1_0 + i * 4,
                                             (const PetscScalar *)buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        