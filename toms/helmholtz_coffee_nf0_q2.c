
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
  static const double IP_J_K_1_3_9[10][10]  = {{0.00238095238095, 0.000396825396825, 0.000396825396825, 0.000396825396825, -0.00238095238095, -0.00238095238095, -0.00238095238095, -0.0015873015873, -0.0015873015873, -0.0015873015873}, 
  {0.000396825396825, 0.00238095238095, 0.000396825396825, 0.000396825396825, -0.00238095238095, -0.0015873015873, -0.0015873015873, -0.00238095238095, -0.00238095238095, -0.0015873015873}, 
  {0.000396825396825, 0.000396825396825, 0.00238095238095, 0.000396825396825, -0.0015873015873, -0.00238095238095, -0.0015873015873, -0.00238095238095, -0.0015873015873, -0.00238095238095}, 
  {0.000396825396825, 0.000396825396825, 0.000396825396825, 0.00238095238095, -0.0015873015873, -0.0015873015873, -0.00238095238095, -0.0015873015873, -0.00238095238095, -0.00238095238095}, 
  {-0.00238095238095, -0.00238095238095, -0.0015873015873, -0.0015873015873, 0.0126984126984, 0.00634920634921, 0.00634920634921, 0.00634920634921, 0.00634920634921, 0.0031746031746}, 
  {-0.00238095238095, -0.0015873015873, -0.00238095238095, -0.0015873015873, 0.00634920634921, 0.0126984126984, 0.00634920634921, 0.00634920634921, 0.0031746031746, 0.00634920634921}, 
  {-0.00238095238095, -0.0015873015873, -0.0015873015873, -0.00238095238095, 0.00634920634921, 0.00634920634921, 0.0126984126984, 0.0031746031746, 0.00634920634921, 0.00634920634921}, 
  {-0.0015873015873, -0.00238095238095, -0.00238095238095, -0.0015873015873, 0.00634920634921, 0.00634920634921, 0.0031746031746, 0.0126984126984, 0.00634920634921, 0.00634920634921}, 
  {-0.0015873015873, -0.00238095238095, -0.0015873015873, -0.00238095238095, 0.00634920634921, 0.0031746031746, 0.00634920634921, 0.00634920634921, 0.0126984126984, 0.00634920634921}, 
  {-0.0015873015873, -0.0015873015873, -0.00238095238095, -0.00238095238095, 0.0031746031746, 0.00634920634921, 0.00634920634921, 0.00634920634921, 0.00634920634921, 0.0126984126984}};
  static const double IP_J_K_1_3_8[10][10]  = {{0.1, 0.0, 0.0333333333333, 0.0, 0.0333333333333, 0.0, 0.0333333333333, -0.0333333333333, -0.133333333333, -0.0333333333333}, 
  {0.0333333333333, 0.0, -0.0333333333333, 0.0, -0.0333333333333, 0.0, 0.1, 0.0333333333333, 0.0, -0.1}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, -0.0333333333333, 0.0, 0.266666666667, 0.0, 0.133333333333, -0.266666666667, -0.0, -0.133333333333}, 
  {0.0333333333333, 0.0, 0.1, 0.0, 0.133333333333, 0.0, 0.133333333333, -0.133333333333, -0.133333333333, -0.133333333333}, 
  {-0.0333333333333, 0.0, 0.0333333333333, 0.0, -0.266666666667, 0.0, -0.133333333333, 0.266666666667, 0.0, 0.133333333333}, 
  {-0.0333333333333, 0.0, -0.1, 0.0, -0.133333333333, 0.0, -0.133333333333, 0.133333333333, 0.133333333333, 0.133333333333}, 
  {-0.133333333333, 0.0, -0.0, 0.0, -0.0, 0.0, -0.133333333333, 0.0, 0.133333333333, 0.133333333333}};
  static const double IP_J_K_1_3_7[10][10]  = {{0.1, 0.0333333333333, 0.0, 0.0, 0.0, 0.0333333333333, 0.0333333333333, -0.0333333333333, -0.0333333333333, -0.133333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, -0.0333333333333, 0.0, 0.0, 0.0, 0.1, -0.0333333333333, -0.1, 0.0333333333333, -0.0}, 
  {0.0333333333333, -0.0333333333333, 0.0, 0.0, 0.0, 0.133333333333, 0.266666666667, -0.133333333333, -0.266666666667, -0.0}, 
  {0.0333333333333, 0.1, 0.0, 0.0, 0.0, 0.133333333333, 0.133333333333, -0.133333333333, -0.133333333333, -0.133333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.133333333333, 0.0, 0.0, 0.0, 0.0, -0.133333333333, -0.0, 0.133333333333, 0.0, 0.133333333333}, 
  {-0.0333333333333, 0.0333333333333, 0.0, 0.0, 0.0, -0.133333333333, -0.266666666667, 0.133333333333, 0.266666666667, 0.0}, 
  {-0.0333333333333, -0.1, 0.0, 0.0, 0.0, -0.133333333333, -0.133333333333, 0.133333333333, 0.133333333333, 0.133333333333}};
  static const double IP_J_K_1_3_6[10][10]  = {{0.1, 0.0, 0.0, 0.0333333333333, 0.0333333333333, 0.0333333333333, 0.0, -0.133333333333, -0.0333333333333, -0.0333333333333}, 
  {0.0333333333333, 0.0, 0.0, -0.0333333333333, -0.0333333333333, 0.1, 0.0, 0.0, 0.0333333333333, -0.1}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, 0.0, 0.1, 0.133333333333, 0.133333333333, 0.0, -0.133333333333, -0.133333333333, -0.133333333333}, 
  {0.0333333333333, 0.0, 0.0, -0.0333333333333, 0.266666666667, 0.133333333333, 0.0, -0.0, -0.266666666667, -0.133333333333}, 
  {-0.0333333333333, 0.0, 0.0, -0.1, -0.133333333333, -0.133333333333, 0.0, 0.133333333333, 0.133333333333, 0.133333333333}, 
  {-0.0333333333333, 0.0, 0.0, 0.0333333333333, -0.266666666667, -0.133333333333, 0.0, 0.0, 0.266666666667, 0.133333333333}, 
  {-0.133333333333, 0.0, 0.0, -0.0, -0.0, -0.133333333333, 0.0, 0.133333333333, 0.0, 0.133333333333}};
  static const double IP_J_K_1_3_5[10][10]  = {{0.1, 0.0, 0.0333333333333, 0.0, 0.0333333333333, 0.0, 0.0333333333333, -0.0333333333333, -0.133333333333, -0.0333333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, -0.0333333333333, 0.0, 0.1, 0.0, -0.0333333333333, -0.1, -0.0, 0.0333333333333}, 
  {0.0333333333333, 0.0, 0.1, 0.0, 0.133333333333, 0.0, 0.133333333333, -0.133333333333, -0.133333333333, -0.133333333333}, 
  {0.0333333333333, 0.0, -0.0333333333333, 0.0, 0.133333333333, 0.0, 0.266666666667, -0.133333333333, -0.0, -0.266666666667}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.133333333333, 0.0, -0.0, 0.0, -0.133333333333, 0.0, 0.0, 0.133333333333, 0.133333333333, 0.0}, 
  {-0.0333333333333, 0.0, -0.1, 0.0, -0.133333333333, 0.0, -0.133333333333, 0.133333333333, 0.133333333333, 0.133333333333}, 
  {-0.0333333333333, 0.0, 0.0333333333333, 0.0, -0.133333333333, 0.0, -0.266666666667, 0.133333333333, 0.0, 0.266666666667}};
  static const double IP_J_K_1_3_4[10][10]  = {{0.1, 0.0333333333333, 0.0, 0.0, 0.0, 0.0333333333333, 0.0333333333333, -0.0333333333333, -0.0333333333333, -0.133333333333}, 
  {0.0333333333333, 0.1, 0.0, 0.0, 0.0, -0.0333333333333, -0.0333333333333, 0.0333333333333, 0.0333333333333, -0.133333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, -0.0333333333333, 0.0, 0.0, 0.0, 0.266666666667, 0.133333333333, -0.266666666667, -0.133333333333, -0.0}, 
  {0.0333333333333, -0.0333333333333, 0.0, 0.0, 0.0, 0.133333333333, 0.266666666667, -0.133333333333, -0.266666666667, -0.0}, 
  {-0.0333333333333, 0.0333333333333, 0.0, 0.0, 0.0, -0.266666666667, -0.133333333333, 0.266666666667, 0.133333333333, 0.0}, 
  {-0.0333333333333, 0.0333333333333, 0.0, 0.0, 0.0, -0.133333333333, -0.266666666667, 0.133333333333, 0.266666666667, 0.0}, 
  {-0.133333333333, -0.133333333333, 0.0, 0.0, 0.0, -0.0, -0.0, 0.0, 0.0, 0.266666666667}};
  static const double IP_J_K_1_3_3[10][10]  = {{0.1, 0.0, 0.0333333333333, 0.0, 0.0333333333333, 0.0, 0.0333333333333, -0.0333333333333, -0.133333333333, -0.0333333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, 0.1, 0.0, -0.0333333333333, 0.0, -0.0333333333333, 0.0333333333333, -0.133333333333, 0.0333333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, -0.0333333333333, 0.0, 0.266666666667, 0.0, 0.133333333333, -0.266666666667, -0.0, -0.133333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, -0.0333333333333, 0.0, 0.133333333333, 0.0, 0.266666666667, -0.133333333333, -0.0, -0.266666666667}, 
  {-0.0333333333333, 0.0, 0.0333333333333, 0.0, -0.266666666667, 0.0, -0.133333333333, 0.266666666667, 0.0, 0.133333333333}, 
  {-0.133333333333, 0.0, -0.133333333333, 0.0, -0.0, 0.0, -0.0, 0.0, 0.266666666667, 0.0}, 
  {-0.0333333333333, 0.0, 0.0333333333333, 0.0, -0.133333333333, 0.0, -0.266666666667, 0.133333333333, 0.0, 0.266666666667}};
  static const double IP_J_K_1_3_2[10][10]  = {{0.1, 0.0, 0.0, 0.0333333333333, 0.0333333333333, 0.0333333333333, 0.0, -0.133333333333, -0.0333333333333, -0.0333333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, 0.0, -0.0333333333333, 0.1, -0.0333333333333, 0.0, -0.0, -0.1, 0.0333333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, 0.0, 0.1, 0.133333333333, 0.133333333333, 0.0, -0.133333333333, -0.133333333333, -0.133333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, 0.0, -0.0333333333333, 0.133333333333, 0.266666666667, 0.0, 0.0, -0.133333333333, -0.266666666667}, 
  {-0.0333333333333, 0.0, 0.0, -0.1, -0.133333333333, -0.133333333333, 0.0, 0.133333333333, 0.133333333333, 0.133333333333}, 
  {-0.133333333333, 0.0, 0.0, -0.0, -0.133333333333, -0.0, 0.0, 0.133333333333, 0.133333333333, 0.0}, 
  {-0.0333333333333, 0.0, 0.0, 0.0333333333333, -0.133333333333, -0.266666666667, 0.0, 0.0, 0.133333333333, 0.266666666667}};
  static const double IP_J_K_1_3_1[10][10]  = {{0.1, 0.0333333333333, 0.0, 0.0, 0.0, 0.0333333333333, 0.0333333333333, -0.0333333333333, -0.0333333333333, -0.133333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, -0.0333333333333, 0.0, 0.0, 0.0, -0.0333333333333, 0.1, 0.0333333333333, -0.1, -0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, -0.0333333333333, 0.0, 0.0, 0.0, 0.266666666667, 0.133333333333, -0.266666666667, -0.133333333333, -0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.1, 0.0, 0.0, 0.0, 0.133333333333, 0.133333333333, -0.133333333333, -0.133333333333, -0.133333333333}, 
  {-0.0333333333333, 0.0333333333333, 0.0, 0.0, 0.0, -0.266666666667, -0.133333333333, 0.266666666667, 0.133333333333, 0.0}, 
  {-0.133333333333, 0.0, 0.0, 0.0, 0.0, -0.0, -0.133333333333, 0.0, 0.133333333333, 0.133333333333}, 
  {-0.0333333333333, -0.1, 0.0, 0.0, 0.0, -0.133333333333, -0.133333333333, 0.133333333333, 0.133333333333, 0.133333333333}};
  static const double IP_J_K_1_3_0[10][10]  = {{0.1, 0.0, 0.0, 0.0333333333333, 0.0333333333333, 0.0333333333333, 0.0, -0.133333333333, -0.0333333333333, -0.0333333333333}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0333333333333, 0.0, 0.0, 0.1, -0.0333333333333, -0.0333333333333, 0.0, -0.133333333333, 0.0333333333333, 0.0333333333333}, 
  {0.0333333333333, 0.0, 0.0, -0.0333333333333, 0.266666666667, 0.133333333333, 0.0, 0.0, -0.266666666667, -0.133333333333}, 
  {0.0333333333333, 0.0, 0.0, -0.0333333333333, 0.133333333333, 0.266666666667, 0.0, -0.0, -0.133333333333, -0.266666666667}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.133333333333, 0.0, 0.0, -0.133333333333, 0.0, -0.0, 0.0, 0.266666666667, -0.0, 0.0}, 
  {-0.0333333333333, 0.0, 0.0, 0.0333333333333, -0.266666666667, -0.133333333333, 0.0, -0.0, 0.266666666667, 0.133333333333}, 
  {-0.0333333333333, 0.0, 0.0, 0.0333333333333, -0.133333333333, -0.266666666667, 0.0, 0.0, 0.133333333333, 0.266666666667}};
  
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
  static const double FE0_D100[14][10]  = {{1.0, -1.0, 0.0, 0.0, 0.0, 2, 2, -2, -2, 0.0}, 
  {1.0, 1.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 0.0, -2}, 
  {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2, 0.0, -2, -2}, 
  {-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 2.0}, 
  {-1.0, -1.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 0.0, 2}, 
  {0.597892939099184, 1.79367881729755, 0.0, 0.0, 0.0, 0.402107060900818, 0.402107060900818, -0.402107060900818, -0.402107060900818, -2.39157175639673}, 
  {-1.79367881729755, -0.597892939099182, 0.0, 0.0, 0.0, 0.402107060900818, 0.402107060900818, -0.402107060900818, -0.402107060900818, 2.39157175639673}, 
  {0.597892939099184, -0.597892939099182, 0.0, 0.0, 0.0, 2.79367881729755, 0.402107060900818, -2.79367881729755, -0.402107060900818, 0.0}, 
  {0.597892939099183, -0.597892939099185, 0.0, 0.0, 0.0, 0.402107060900818, 2.79367881729755, -0.402107060900818, -2.79367881729755, 0.0}, 
  {-0.257491493972768, -0.772474481918307, 0.0, 0.0, 0.0, 1.25749149397277, 1.25749149397277, -1.25749149397277, -1.25749149397277, 1.02996597589107}, 
  {0.772474481918307, 0.257491493972768, 0.0, 0.0, 0.0, 1.25749149397277, 1.25749149397277, -1.25749149397277, -1.25749149397277, -1.02996597589108}, 
  {-0.257491493972768, 0.257491493972768, 0.0, 0.0, 0.0, 0.227525518081694, 1.25749149397277, -0.227525518081694, -1.25749149397277, 0.0}, 
  {-0.257491493972769, 0.257491493972768, 0.0, 0.0, 0.0, 1.25749149397277, 0.227525518081694, -1.25749149397277, -0.227525518081694, 0.0}};
  static const double FE0_D001[14][10]  = {{1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, -2, -2.0, 0.0}, 
  {1.0, 0.0, 0.0, 1.0, 0.0, 2.00000000000001, 0.0, -2, 0.0, -2.0}, 
  {1.0, 0.0, 0.0, -1.0, 2, 2.00000000000001, 0.0, 0.0, -2.0, -2}, 
  {-1.0, 0.0, 0.0, -1.0, 0.0, 2.00000000000001, 0.0, 2.0, 0.0, -2}, 
  {-1.0, 0.0, 0.0, -1.0, 2, 0.0, 0.0, 2, -2, 0.0}, 
  {-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099184, 0.0, 0.0, -0.597892939099186, 0.40210706090082, 2.79367881729756, 0.0, 0.0, -0.402107060900818, -2.79367881729755}, 
  {-1.79367881729755, 0.0, 0.0, -0.597892939099187, 0.40210706090082, 0.402107060900823, 0.0, 2.39157175639673, -0.402107060900818, -0.402107060900818}, 
  {0.597892939099184, 0.0, 0.0, 1.79367881729754, 0.402107060900815, 0.402107060900826, 0.0, -2.39157175639673, -0.402107060900816, -0.402107060900819}, 
  {0.597892939099183, 0.0, 0.0, -0.597892939099186, 2.79367881729755, 0.402107060900822, 0.0, 0.0, -2.79367881729755, -0.402107060900818}, 
  {-0.257491493972768, 0.0, 0.0, 0.257491493972765, 1.25749149397277, 0.227525518081701, 0.0, 0.0, -1.25749149397277, -0.227525518081694}, 
  {0.772474481918307, 0.0, 0.0, 0.257491493972765, 1.25749149397277, 1.25749149397278, 0.0, -1.02996597589108, -1.25749149397277, -1.25749149397277}, 
  {-0.257491493972768, 0.0, 0.0, -0.772474481918311, 1.25749149397277, 1.25749149397278, 0.0, 1.02996597589107, -1.25749149397277, -1.25749149397277}, 
  {-0.257491493972768, 0.0, 0.0, 0.257491493972765, 0.227525518081694, 1.25749149397278, 0.0, 0.0, -0.227525518081694, -1.25749149397277}};
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
  static const double FE0_D010[14][10]  = {{1.0, 0.0, 1.0, 0.0, 2, 0.0, 0.0, -2, -2, 0.0}, 
  {1.0, 0.0, -1.0, 0.0, 2, 0.0, 2, -2, 0.0, -2.0}, 
  {1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.00000000000001, 0.0, -2.0, -2.0}, 
  {-1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 2.00000000000001, 0.0, 2, -2}, 
  {-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, 0.0, -1.0, 0.0, 2, 0.0, 0.0, -2, 2.0, 0.0}, 
  {0.597892939099184, 0.0, -0.597892939099185, 0.0, 0.402107060900819, 0.0, 2.79367881729755, -0.402107060900818, 0.0, -2.79367881729755}, 
  {-1.79367881729755, 0.0, -0.597892939099185, 0.0, 0.402107060900819, 0.0, 0.402107060900823, -0.402107060900818, 2.39157175639673, -0.402107060900818}, 
  {0.597892939099185, 0.0, -0.597892939099184, 0.0, 2.79367881729755, 0.0, 0.402107060900821, -2.79367881729755, 0.0, -0.402107060900818}, 
  {0.597892939099184, 0.0, 1.79367881729754, 0.0, 0.402107060900819, 0.0, 0.402107060900825, -0.402107060900818, -2.39157175639673, -0.402107060900819}, 
  {-0.257491493972768, 0.0, 0.257491493972766, 0.0, 1.25749149397277, 0.0, 0.227525518081699, -1.25749149397277, 0.0, -0.227525518081694}, 
  {0.772474481918307, 0.0, 0.257491493972766, 0.0, 1.25749149397277, 0.0, 1.25749149397277, -1.25749149397277, -1.02996597589108, -1.25749149397277}, 
  {-0.257491493972767, 0.0, 0.257491493972765, 0.0, 0.227525518081694, 0.0, 1.25749149397278, -0.227525518081694, 0.0, -1.25749149397277}, 
  {-0.257491493972767, 0.0, -0.772474481918308, 0.0, 1.25749149397277, 0.0, 1.25749149397277, -1.25749149397277, 1.02996597589107, -1.25749149397277}};
  double CONST_1_1_0  = ((K[8] * K[5] * det) + (K[7] * K[4] * det) + (K[6] * K[3] * det));
  double CONST_1_1_1  = ((K[2] * K[8] * det) + (K[1] * K[7] * det) + (K[0] * K[6] * det));
  double CONST_1_1_2  = ((K[2] * K[2] * det) + (K[1] * K[1] * det) + (K[0] * K[0] * det));
  double CONST_1_1_3  = ((K[5] * K[2] * det) + (K[4] * K[1] * det) + (K[3] * K[0] * det));
  double CONST_1_1_4  = ((K[2] * K[5] * det) + (K[1] * K[4] * det) + (K[0] * K[3] * det));
  double CONST_1_1_5  = ((K[8] * K[2] * det) + (K[7] * K[1] * det) + (K[6] * K[0] * det));
  double CONST_1_1_6  = ((K[5] * K[5] * det) + (K[4] * K[4] * det) + (K[3] * K[3] * det));
  double CONST_1_1_7  = ((K[8] * K[8] * det) + (K[7] * K[7] * det) + (K[6] * K[6] * det));
  double CONST_1_1_8  = ((K[5] * K[8] * det) + (K[4] * K[7] * det) + (K[3] * K[6] * det));
  
  
  for (int j  = 0; j < 10; j += 1)
  {
    
    for (int k  = 0; k < 10; k += 1)
    {
      A[j][k] += (((IP_J_K_1_3_9[j][k] * det) + (IP_J_K_1_3_4[j][k] * CONST_1_1_2) + (IP_J_K_1_3_1[j][k] * CONST_1_1_4) + (IP_J_K_1_3_7[j][k] * CONST_1_1_1) + (IP_J_K_1_3_8[j][k] * CONST_1_1_3) + (IP_J_K_1_3_3[j][k] * CONST_1_1_6) + (IP_J_K_1_3_5[j][k] * CONST_1_1_8) + (IP_J_K_1_3_6[j][k] * CONST_1_1_5) + (IP_J_K_1_3_2[j][k] * CONST_1_1_0) + (IP_J_K_1_3_0[j][k] * CONST_1_1_7)));
      
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
        
        