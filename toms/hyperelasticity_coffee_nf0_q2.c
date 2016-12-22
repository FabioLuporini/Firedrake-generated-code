
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[30][30] , double** coordinate_dofs , double** w0 , double* c1 , double* c2 )
{
  // This code is generated visiting a COFFEE AST
  
  double (*w1)[1] = (double (*)[1])c1;
  
  double (*w2)[1] = (double (*)[1])c2;
  
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
  static const double FE1_D010[14][10]  = {{1.0, 0.0, 1.0, 0.0, 2, 0.0, 0.0, -2, -2, 0.0}, 
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
  static const double FE1_D001[14][10]  = {{1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, -2, -2.0, 0.0}, 
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
  static const double FE1_D100[14][10]  = {{1.0, -1.0, 0.0, 0.0, 0.0, 2, 2, -2, -2, 0.0}, 
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
  static const double FE2[14][1]  = {{1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}, 
  {1.0}};
  double F9  = 0.0;
  double F10  = 0.0;
  
  
  for (int r  = 0; r < 1; r += 1)
  {
    F9 += (w1[0][0] * FE2[0][0]);
    F10 += (w2[0][0] * FE2[0][0]);
    
  }
  
  
  for (int ip  = 0; ip < 14; ip += 1)
  {
    double F0  = 0.0;
    double F1  = 0.0;
    double F2  = 0.0;
    double F3  = 0.0;
    double F4  = 0.0;
    double F5  = 0.0;
    double F6  = 0.0;
    double F7  = 0.0;
    double F8  = 0.0;
    
    for (int r  = 0; r < 10; r += 1)
    {
      F5 += (w0[r+10][0] * FE1_D001[ip][r]);
      F4 += (w0[r+10][0] * FE1_D010[ip][r]);
      F3 += (w0[r+10][0] * FE1_D100[ip][r]);
      F8 += (w0[r+20][0] * FE1_D001[ip][r]);
      F7 += (w0[r+20][0] * FE1_D010[ip][r]);
      F6 += (w0[r+20][0] * FE1_D100[ip][r]);
      F2 += (w0[r][0] * FE1_D001[ip][r]);
      F1 += (w0[r][0] * FE1_D010[ip][r]);
      F0 += (w0[r][0] * FE1_D100[ip][r]);
      
    }
    double IP_0_1_0  = (2 * 0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0))));
    double IP_0_1_1  = (2 * 0.5 * ((((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + (((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0))));
    double IP_0_1_2  = ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0);
    double IP_0_1_3  = ((F9 * 2 * 0.5 * ((((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0)) + -1.0)) + (F10 * ((0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[8] * F5) + (K[5] * F4) + (K[2] * F3))) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[8] * F2) + (K[5] * F1) + (K[2] * F0))) + (((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + -1.0)) + (0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[7] * F8) + (K[4] * F7) + (K[1] * F6))) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + -1.0)) + (0.5 * ((((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0)) + -1.0)))));
    double IP_0_1_4  = ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0);
    double IP_0_1_5  = ((K[6] * F5) + (K[3] * F4) + (K[0] * F3));
    double IP_0_1_6  = ((F9 * 2 * 0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[7] * F8) + (K[4] * F7) + (K[1] * F6))) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + -1.0)) + (F10 * ((0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[8] * F5) + (K[5] * F4) + (K[2] * F3))) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[8] * F2) + (K[5] * F1) + (K[2] * F0))) + (((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + -1.0)) + (0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[7] * F8) + (K[4] * F7) + (K[1] * F6))) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + -1.0)) + (0.5 * ((((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0)) + -1.0)))));
    double IP_0_1_7  = ((K[6] * F8) + (K[3] * F7) + (K[0] * F6));
    double IP_0_1_8  = ((K[8] * F2) + (K[5] * F1) + (K[2] * F0));
    double IP_0_1_9  = ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0);
    double IP_0_1_10  = (2 * 0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0))));
    double IP_0_1_11  = ((K[8] * F5) + (K[5] * F4) + (K[2] * F3));
    double IP_0_1_12  = ((K[7] * F8) + (K[4] * F7) + (K[1] * F6));
    double IP_0_1_13  = ((K[7] * F2) + (K[4] * F1) + (K[1] * F0));
    double IP_0_1_14  = ((F9 * 2 * 0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[8] * F5) + (K[5] * F4) + (K[2] * F3))) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[8] * F2) + (K[5] * F1) + (K[2] * F0))) + (((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + -1.0)) + (F10 * ((0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[8] * F5) + (K[5] * F4) + (K[2] * F3))) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[8] * F2) + (K[5] * F1) + (K[2] * F0))) + (((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + -1.0)) + (0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[7] * F8) + (K[4] * F7) + (K[1] * F6))) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + -1.0)) + (0.5 * ((((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0)) + -1.0)))));
    double IP_J_0_3_0[30]  = {0.0};
    double IP_J_0_3_1[30]  = {0.0};
    double IP_J_0_3_2[30]  = {0.0};
    double IP_J_0_3_3[30]  = {0.0};
    double IP_J_0_3_4[30]  = {0.0};
    double IP_J_0_3_5[30]  = {0.0};
    double IP_J_0_3_6[30]  = {0.0};
    double IP_J_0_3_7[30]  = {0.0};
    double IP_J_0_3_8[30]  = {0.0};
    double IP_K_0_3_0_EXP_0_2[30]  = {0.0};
    double IP_K_0_3_1_EXP_0_4[30]  = {0.0};
    double IP_K_0_3_2_EXP_0_5[30]  = {0.0};
    double IP_K_0_3_3_EXP_0_7[30]  = {0.0};
    double IP_K_0_3_4_EXP_0_8[30]  = {0.0};
    double IP_K_0_3_5_EXP_0_6[30]  = {0.0};
    double IP_K_0_3_6_EXP_0_1[30]  = {0.0};
    double IP_K_0_3_7_EXP_0_3[30]  = {0.0};
    double IP_K_0_3_8_EXP_0_0[30]  = {0.0};
    double IP_0_4_0  = ((IP_0_1_13 * IP_0_1_4 * F9) + (IP_0_1_8 * F10 * IP_0_1_12));
    double IP_0_4_1  = ((IP_0_1_2 * IP_0_1_12 * F9) + (IP_0_1_13 * F10 * IP_0_1_7));
    double IP_0_4_2  = ((IP_0_1_2 * IP_0_1_9 * F9) + (IP_0_1_13 * F10 * IP_0_1_5));
    double IP_0_4_3  = ((IP_0_1_9 * IP_0_1_7 * F9) + (IP_0_1_5 * F10 * IP_0_1_12));
    double IP_0_4_4  = ((IP_0_1_13 * IP_0_1_13 * F9) + (IP_0_1_2 * IP_0_1_2 * F9) + (IP_0_1_8 * IP_0_1_8 * ((2.0 * F9) + F10)) + IP_0_1_14);
    double IP_0_4_5  = ((F9 * (IP_0_1_0 + (IP_0_1_13 * IP_0_1_2))) + (IP_0_1_2 * F10 * IP_0_1_13));
    double IP_0_4_6  = ((IP_0_1_12 * IP_0_1_12 * F9) + (IP_0_1_7 * IP_0_1_7 * F9) + IP_0_1_14 + (IP_0_1_4 * IP_0_1_4 * ((2.0 * F9) + F10)));
    double IP_0_4_7  = ((F9 * ((IP_0_1_2 * IP_0_1_5) + (IP_0_1_13 * IP_0_1_9) + (2.0 * IP_0_1_8 * IP_0_1_11))) + (IP_0_1_8 * F10 * IP_0_1_11));
    double IP_0_4_8  = ((IP_0_1_8 * IP_0_1_7 * F9) + (IP_0_1_2 * F10 * IP_0_1_4));
    double IP_0_4_9  = ((IP_0_1_12 * IP_0_1_8 * F9) + (IP_0_1_4 * F10 * IP_0_1_13));
    double IP_0_4_10  = ((F9 * ((IP_0_1_13 * IP_0_1_12) + (IP_0_1_2 * IP_0_1_7) + (2.0 * IP_0_1_8 * IP_0_1_4))) + (IP_0_1_8 * F10 * IP_0_1_4));
    double IP_0_4_11  = ((F9 * ((IP_0_1_2 * IP_0_1_7) + (IP_0_1_8 * IP_0_1_4) + (2.0 * IP_0_1_13 * IP_0_1_12))) + (IP_0_1_13 * F10 * IP_0_1_12));
    double IP_0_4_12  = ((F9 * ((IP_0_1_8 * IP_0_1_11) + (IP_0_1_13 * IP_0_1_9) + (2.0 * IP_0_1_2 * IP_0_1_5))) + (IP_0_1_2 * F10 * IP_0_1_5));
    double IP_0_4_13  = ((F9 * ((IP_0_1_9 * IP_0_1_12) + (IP_0_1_11 * IP_0_1_4) + (2.0 * IP_0_1_5 * IP_0_1_7))) + (IP_0_1_5 * F10 * IP_0_1_7));
    double IP_0_4_14  = ((IP_0_1_13 * IP_0_1_11 * F9) + (IP_0_1_8 * F10 * IP_0_1_9));
    double IP_0_4_15  = ((F9 * ((IP_0_1_11 * IP_0_1_8) + (IP_0_1_5 * IP_0_1_2) + (2.0 * IP_0_1_9 * IP_0_1_13))) + (IP_0_1_9 * F10 * IP_0_1_13));
    double IP_0_4_16  = ((IP_0_1_7 * IP_0_1_12 * F9) + (IP_0_1_0 * F9) + (IP_0_1_12 * IP_0_1_7 * F10));
    double IP_0_4_17  = ((F9 * ((IP_0_1_4 * IP_0_1_11) + (IP_0_1_12 * IP_0_1_9) + (2.0 * IP_0_1_7 * IP_0_1_5))) + (IP_0_1_7 * F10 * IP_0_1_5));
    double IP_0_4_18  = ((IP_0_1_4 * IP_0_1_9 * F9) + (IP_0_1_12 * F10 * IP_0_1_11));
    double IP_0_4_19  = ((IP_0_1_2 * IP_0_1_13 * F9) + (IP_0_1_0 * F9) + (IP_0_1_13 * IP_0_1_2 * F10));
    double IP_0_4_20  = ((IP_0_1_11 * IP_0_1_5 * F9) + (IP_0_1_10 * F9) + (IP_0_1_5 * IP_0_1_11 * F10));
    double IP_0_4_21  = ((IP_0_1_4 * IP_0_1_7 * F9) + (IP_0_1_10 * F9) + (IP_0_1_7 * IP_0_1_4 * F10));
    double IP_0_4_22  = ((IP_0_1_7 * IP_0_1_11 * F9) + (IP_0_1_4 * F10 * IP_0_1_5));
    double IP_0_4_23  = ((IP_0_1_9 * IP_0_1_11 * F9) + (IP_0_1_1 * F9) + (IP_0_1_11 * IP_0_1_9 * F10));
    double IP_0_4_24  = ((IP_0_1_11 * IP_0_1_12 * F9) + (IP_0_1_9 * F10 * IP_0_1_4));
    double IP_0_4_25  = ((IP_0_1_8 * IP_0_1_12 * F9) + (IP_0_1_13 * F10 * IP_0_1_4));
    double IP_0_4_26  = ((IP_0_1_5 * IP_0_1_12 * F9) + (IP_0_1_9 * F10 * IP_0_1_7));
    double IP_0_4_27  = ((IP_0_1_7 * IP_0_1_13 * F9) + (IP_0_1_12 * F10 * IP_0_1_2));
    double IP_0_4_28  = ((IP_0_1_7 * IP_0_1_9 * F9) + (IP_0_1_12 * F10 * IP_0_1_5));
    double IP_0_4_29  = ((F9 * ((IP_0_1_5 * IP_0_1_7) + (IP_0_1_11 * IP_0_1_4) + (2.0 * IP_0_1_9 * IP_0_1_12))) + (IP_0_1_9 * F10 * IP_0_1_12));
    double IP_0_4_30  = ((IP_0_1_9 * IP_0_1_2 * F9) + (IP_0_1_5 * F10 * IP_0_1_13));
    double IP_0_4_31  = ((IP_0_1_5 * IP_0_1_13 * F9) + (IP_0_1_9 * F10 * IP_0_1_2));
    double IP_0_4_32  = ((IP_0_1_12 * IP_0_1_12 * F9) + (IP_0_1_4 * IP_0_1_4 * F9) + (IP_0_1_7 * IP_0_1_7 * ((2.0 * F9) + F10)) + IP_0_1_3);
    double IP_0_4_33  = ((IP_0_1_13 * IP_0_1_5 * F9) + (IP_0_1_2 * F10 * IP_0_1_9));
    double IP_0_4_34  = ((IP_0_1_7 * IP_0_1_7 * F9) + (IP_0_1_4 * IP_0_1_4 * F9) + (IP_0_1_12 * IP_0_1_12 * ((2.0 * F9) + F10)) + IP_0_1_6);
    double IP_0_4_35  = ((IP_0_1_13 * IP_0_1_8 * F9) + (IP_0_1_1 * F9) + (IP_0_1_8 * IP_0_1_13 * F10));
    double IP_0_4_36  = ((IP_0_1_11 * IP_0_1_7 * F9) + (IP_0_1_5 * F10 * IP_0_1_4));
    double IP_0_4_37  = ((IP_0_1_4 * IP_0_1_12 * F9) + (IP_0_1_1 * F9) + (IP_0_1_12 * IP_0_1_4 * F10));
    double IP_0_4_38  = ((F9 * ((IP_0_1_4 * IP_0_1_8) + (IP_0_1_12 * IP_0_1_13) + (2.0 * IP_0_1_7 * IP_0_1_2))) + (IP_0_1_7 * F10 * IP_0_1_2));
    double IP_0_4_39  = ((IP_0_1_8 * IP_0_1_13 * F9) + (IP_0_1_1 * F9) + (IP_0_1_13 * IP_0_1_8 * F10));
    double IP_0_4_40  = ((F9 * ((IP_0_1_13 * IP_0_1_12) + (IP_0_1_8 * IP_0_1_4) + (2.0 * IP_0_1_2 * IP_0_1_7))) + (IP_0_1_2 * F10 * IP_0_1_7));
    double IP_0_4_41  = ((F9 * ((IP_0_1_9 * IP_0_1_12) + (IP_0_1_5 * IP_0_1_7) + (2.0 * IP_0_1_11 * IP_0_1_4))) + (IP_0_1_11 * F10 * IP_0_1_4));
    double IP_0_4_42  = ((IP_0_1_9 * IP_0_1_8 * F9) + (IP_0_1_11 * F10 * IP_0_1_13));
    double IP_0_4_43  = ((IP_0_1_12 * IP_0_1_7 * F9) + (IP_0_1_0 * F9) + (IP_0_1_7 * IP_0_1_12 * F10));
    double IP_0_4_44  = ((F9 * ((IP_0_1_11 * IP_0_1_8) + (IP_0_1_9 * IP_0_1_13) + (2.0 * IP_0_1_5 * IP_0_1_2))) + (IP_0_1_5 * F10 * IP_0_1_2));
    double IP_0_4_45  = ((IP_0_1_8 * IP_0_1_8 * F9) + (IP_0_1_2 * IP_0_1_2 * F9) + (IP_0_1_13 * IP_0_1_13 * ((2.0 * F9) + F10)) + IP_0_1_6);
    double IP_0_4_46  = ((IP_0_1_5 * IP_0_1_5 * F9) + (IP_0_1_9 * IP_0_1_9 * F9) + (IP_0_1_11 * IP_0_1_11 * ((2.0 * F9) + F10)) + IP_0_1_14);
    double IP_0_4_47  = ((IP_0_1_11 * IP_0_1_11 * F9) + (IP_0_1_9 * IP_0_1_9 * F9) + (IP_0_1_5 * IP_0_1_5 * ((2.0 * F9) + F10)) + IP_0_1_3);
    double IP_0_4_48  = ((IP_0_1_5 * IP_0_1_4 * F9) + (IP_0_1_11 * F10 * IP_0_1_7));
    double IP_0_4_49  = ((F9 * (IP_0_1_10 + (IP_0_1_8 * IP_0_1_2))) + (IP_0_1_2 * F10 * IP_0_1_8));
    double IP_0_4_50  = ((IP_0_1_8 * IP_0_1_8 * F9) + (IP_0_1_13 * IP_0_1_13 * F9) + IP_0_1_3 + (IP_0_1_2 * IP_0_1_2 * ((2.0 * F9) + F10)));
    double IP_0_4_51  = ((IP_0_1_9 * IP_0_1_4 * F9) + (IP_0_1_11 * F10 * IP_0_1_12));
    double IP_0_4_52  = ((IP_0_1_11 * IP_0_1_13 * F9) + (IP_0_1_9 * F10 * IP_0_1_8));
    double IP_0_4_53  = ((IP_0_1_12 * IP_0_1_2 * F9) + (IP_0_1_7 * F10 * IP_0_1_13));
    double IP_0_4_54  = ((F9 * ((IP_0_1_4 * IP_0_1_11) + (IP_0_1_7 * IP_0_1_5) + (2.0 * IP_0_1_12 * IP_0_1_9))) + (IP_0_1_12 * F10 * IP_0_1_9));
    double IP_0_4_55  = ((IP_0_1_11 * IP_0_1_11 * F9) + (IP_0_1_5 * IP_0_1_5 * F9) + IP_0_1_6 + (IP_0_1_9 * IP_0_1_9 * ((2.0 * F9) + F10)));
    double IP_0_4_56  = ((IP_0_1_8 * IP_0_1_9 * F9) + (IP_0_1_13 * F10 * IP_0_1_11));
    double IP_0_4_57  = ((IP_0_1_8 * IP_0_1_5 * F9) + (IP_0_1_2 * F10 * IP_0_1_11));
    double IP_0_4_58  = ((IP_0_1_2 * IP_0_1_4 * F9) + (IP_0_1_8 * F10 * IP_0_1_7));
    double IP_0_4_59  = ((IP_0_1_11 * IP_0_1_2 * F9) + (IP_0_1_5 * F10 * IP_0_1_8));
    double IP_0_4_60  = ((IP_0_1_5 * IP_0_1_8 * F9) + (IP_0_1_11 * F10 * IP_0_1_2));
    double IP_0_4_61  = ((IP_0_1_12 * IP_0_1_5 * F9) + (IP_0_1_7 * F10 * IP_0_1_9));
    double IP_0_4_62  = ((F9 * (IP_0_1_10 + (IP_0_1_7 * IP_0_1_4))) + (IP_0_1_4 * F10 * IP_0_1_7));
    double IP_0_4_63  = ((F9 * (IP_0_1_1 + (IP_0_1_12 * IP_0_1_4))) + (IP_0_1_4 * F10 * IP_0_1_12));
    double IP_0_4_64  = ((F9 * ((IP_0_1_9 * IP_0_1_13) + (IP_0_1_5 * IP_0_1_2) + (2.0 * IP_0_1_11 * IP_0_1_8))) + (IP_0_1_11 * F10 * IP_0_1_8));
    double IP_0_4_65  = ((F9 * (IP_0_1_0 + (IP_0_1_5 * IP_0_1_9))) + (IP_0_1_9 * F10 * IP_0_1_5));
    double IP_0_4_66  = ((IP_0_1_5 * IP_0_1_11 * F9) + (IP_0_1_10 * F9) + (IP_0_1_11 * IP_0_1_5 * F10));
    double IP_0_4_67  = ((IP_0_1_4 * IP_0_1_5 * F9) + (IP_0_1_7 * F10 * IP_0_1_11));
    double IP_0_4_68  = ((F9 * ((IP_0_1_12 * IP_0_1_13) + (IP_0_1_7 * IP_0_1_2) + (2.0 * IP_0_1_4 * IP_0_1_8))) + (IP_0_1_4 * F10 * IP_0_1_8));
    double IP_0_4_69  = ((IP_0_1_9 * IP_0_1_5 * F9) + (IP_0_1_0 * F9) + (IP_0_1_5 * IP_0_1_9 * F10));
    double IP_0_4_70  = ((F9 * (IP_0_1_1 + (IP_0_1_11 * IP_0_1_9))) + (IP_0_1_9 * F10 * IP_0_1_11));
    double IP_0_4_71  = ((IP_0_1_13 * IP_0_1_7 * F9) + (IP_0_1_2 * F10 * IP_0_1_12));
    double IP_0_4_72  = ((IP_0_1_7 * IP_0_1_8 * F9) + (IP_0_1_4 * F10 * IP_0_1_2));
    double IP_0_4_73  = ((F9 * ((IP_0_1_7 * IP_0_1_5) + (IP_0_1_12 * IP_0_1_9) + (2.0 * IP_0_1_4 * IP_0_1_11))) + (IP_0_1_4 * F10 * IP_0_1_11));
    double IP_0_4_74  = ((IP_0_1_2 * IP_0_1_8 * F9) + (IP_0_1_10 * F9) + (IP_0_1_8 * IP_0_1_2 * F10));
    double IP_0_4_75  = ((IP_0_1_4 * IP_0_1_13 * F9) + (IP_0_1_12 * F10 * IP_0_1_8));
    double IP_0_4_76  = ((F9 * ((IP_0_1_8 * IP_0_1_11) + (IP_0_1_2 * IP_0_1_5) + (2.0 * IP_0_1_13 * IP_0_1_9))) + (IP_0_1_13 * F10 * IP_0_1_9));
    double IP_0_4_77  = ((IP_0_1_4 * IP_0_1_2 * F9) + (IP_0_1_7 * F10 * IP_0_1_8));
    double IP_0_4_78  = ((IP_0_1_12 * IP_0_1_11 * F9) + (IP_0_1_4 * F10 * IP_0_1_9));
    double IP_0_4_79  = ((F9 * ((IP_0_1_4 * IP_0_1_8) + (IP_0_1_7 * IP_0_1_2) + (2.0 * IP_0_1_12 * IP_0_1_13))) + (IP_0_1_12 * F10 * IP_0_1_13));
    double IP_0_4_80  = ((IP_0_1_2 * IP_0_1_11 * F9) + (IP_0_1_8 * F10 * IP_0_1_5));
    double IP_J_0_6_0[30]  = {0.0};
    double IP_J_0_6_1[30]  = {0.0};
    double IP_J_0_6_2[30]  = {0.0};
    double IP_J_0_6_3[30]  = {0.0};
    double IP_J_0_6_4[30]  = {0.0};
    double IP_J_0_6_5[30]  = {0.0};
    double IP_J_0_6_6[30]  = {0.0};
    double IP_J_0_6_7[30]  = {0.0};
    double IP_J_0_6_8[30]  = {0.0};
    
    for (int j  = 0; j < 10; j += 1)
    {
      IP_J_0_3_0[j+10] = (K[0] * FE1_D100[ip][j]) + (K[3] * FE1_D010[ip][j]) + (K[6] * FE1_D001[ip][j]);
      IP_J_0_3_1[j+20] = (K[2] * FE1_D100[ip][j]) + (K[5] * FE1_D010[ip][j]) + (K[8] * FE1_D001[ip][j]);
      IP_J_0_3_2[j] = (K[2] * FE1_D100[ip][j]) + (K[5] * FE1_D010[ip][j]) + (K[8] * FE1_D001[ip][j]);
      IP_J_0_3_3[j+10] = (K[2] * FE1_D100[ip][j]) + (K[5] * FE1_D010[ip][j]) + (K[8] * FE1_D001[ip][j]);
      IP_J_0_3_4[j+10] = (K[1] * FE1_D100[ip][j]) + (K[4] * FE1_D010[ip][j]) + (K[7] * FE1_D001[ip][j]);
      IP_J_0_3_5[j] = (K[1] * FE1_D100[ip][j]) + (K[4] * FE1_D010[ip][j]) + (K[7] * FE1_D001[ip][j]);
      IP_J_0_3_6[j+20] = (K[0] * FE1_D100[ip][j]) + (K[3] * FE1_D010[ip][j]) + (K[6] * FE1_D001[ip][j]);
      IP_J_0_3_7[j] = (K[0] * FE1_D100[ip][j]) + (K[3] * FE1_D010[ip][j]) + (K[6] * FE1_D001[ip][j]);
      IP_J_0_3_8[j+20] = (K[1] * FE1_D100[ip][j]) + (K[4] * FE1_D010[ip][j]) + (K[7] * FE1_D001[ip][j]);
      IP_K_0_3_8_EXP_0_0[j+20] = ((FE1_D100[ip][j] * K[2]) + (FE1_D010[ip][j] * K[5]) + (FE1_D001[ip][j] * K[8])) * det * W14[ip];
      IP_K_0_3_6_EXP_0_1[j] = ((FE1_D100[ip][j] * K[2]) + (FE1_D010[ip][j] * K[5]) + (FE1_D001[ip][j] * K[8])) * det * W14[ip];
      IP_K_0_3_0_EXP_0_2[j+10] = ((FE1_D100[ip][j] * K[0]) + (FE1_D010[ip][j] * K[3]) + (FE1_D001[ip][j] * K[6])) * det * W14[ip];
      IP_K_0_3_7_EXP_0_3[j] = ((FE1_D100[ip][j] * K[1]) + (FE1_D010[ip][j] * K[4]) + (FE1_D001[ip][j] * K[7])) * det * W14[ip];
      IP_K_0_3_1_EXP_0_4[j+10] = ((FE1_D100[ip][j] * K[1]) + (FE1_D010[ip][j] * K[4]) + (FE1_D001[ip][j] * K[7])) * det * W14[ip];
      IP_K_0_3_2_EXP_0_5[j+20] = ((FE1_D100[ip][j] * K[1]) + (FE1_D010[ip][j] * K[4]) + (FE1_D001[ip][j] * K[7])) * det * W14[ip];
      IP_K_0_3_5_EXP_0_6[j] = ((FE1_D100[ip][j] * K[0]) + (FE1_D010[ip][j] * K[3]) + (FE1_D001[ip][j] * K[6])) * det * W14[ip];
      IP_K_0_3_3_EXP_0_7[j+20] = ((FE1_D100[ip][j] * K[0]) + (FE1_D010[ip][j] * K[3]) + (FE1_D001[ip][j] * K[6])) * det * W14[ip];
      IP_K_0_3_4_EXP_0_8[j+10] = ((FE1_D100[ip][j] * K[2]) + (FE1_D010[ip][j] * K[5]) + (FE1_D001[ip][j] * K[8])) * det * W14[ip];
      
    }
    
    for (int j  = 0; j < 30; j += 1)
    {
      IP_J_0_6_0[j] = (IP_J_0_3_2[j] * IP_0_4_64) + (IP_J_0_3_7[j] * IP_0_4_60) + (IP_J_0_3_3[j] * IP_0_4_46) + (IP_J_0_3_6[j] * IP_0_4_48) + (IP_J_0_3_8[j] * IP_0_4_51) + (IP_J_0_3_4[j] * IP_0_4_23) + (IP_J_0_3_5[j] * IP_0_4_42) + (IP_J_0_3_0[j] * IP_0_4_66) + (IP_J_0_3_1[j] * IP_0_4_41);
      IP_J_0_6_1[j] = (IP_J_0_3_2[j] * IP_0_4_68) + (IP_J_0_3_7[j] * IP_0_4_72) + (IP_J_0_3_3[j] * IP_0_4_73) + (IP_J_0_3_6[j] * IP_0_4_62) + (IP_J_0_3_8[j] * IP_0_4_63) + (IP_J_0_3_4[j] * IP_0_4_78) + (IP_J_0_3_5[j] * IP_0_4_9) + (IP_J_0_3_0[j] * IP_0_4_22) + (IP_J_0_3_1[j] * IP_0_4_6);
      IP_J_0_6_2[j] = (IP_J_0_3_2[j] * IP_0_4_39) + (IP_J_0_3_7[j] * IP_0_4_19) + (IP_J_0_3_3[j] * IP_0_4_56) + (IP_J_0_3_6[j] * IP_0_4_1) + (IP_J_0_3_8[j] * IP_0_4_11) + (IP_J_0_3_4[j] * IP_0_4_76) + (IP_J_0_3_5[j] * IP_0_4_45) + (IP_J_0_3_0[j] * IP_0_4_2) + (IP_J_0_3_1[j] * IP_0_4_25);
      IP_J_0_6_3[j] = (IP_J_0_3_2[j] * IP_0_4_4) + (IP_J_0_3_7[j] * IP_0_4_74) + (IP_J_0_3_3[j] * IP_0_4_7) + (IP_J_0_3_6[j] * IP_0_4_58) + (IP_J_0_3_8[j] * IP_0_4_0) + (IP_J_0_3_4[j] * IP_0_4_14) + (IP_J_0_3_5[j] * IP_0_4_35) + (IP_J_0_3_0[j] * IP_0_4_80) + (IP_J_0_3_1[j] * IP_0_4_10);
      IP_J_0_6_4[j] = (IP_J_0_3_2[j] * IP_0_4_75) + (IP_J_0_3_7[j] * IP_0_4_27) + (IP_J_0_3_3[j] * IP_0_4_18) + (IP_J_0_3_6[j] * IP_0_4_16) + (IP_J_0_3_8[j] * IP_0_4_34) + (IP_J_0_3_4[j] * IP_0_4_54) + (IP_J_0_3_5[j] * IP_0_4_79) + (IP_J_0_3_0[j] * IP_0_4_28) + (IP_J_0_3_1[j] * IP_0_4_37);
      IP_J_0_6_5[j] = (IP_J_0_3_2[j] * IP_0_4_52) + (IP_J_0_3_7[j] * IP_0_4_31) + (IP_J_0_3_3[j] * IP_0_4_70) + (IP_J_0_3_6[j] * IP_0_4_26) + (IP_J_0_3_8[j] * IP_0_4_29) + (IP_J_0_3_4[j] * IP_0_4_55) + (IP_J_0_3_5[j] * IP_0_4_15) + (IP_J_0_3_0[j] * IP_0_4_65) + (IP_J_0_3_1[j] * IP_0_4_24);
      IP_J_0_6_6[j] = (IP_J_0_3_2[j] * IP_0_4_49) + (IP_J_0_3_7[j] * IP_0_4_50) + (IP_J_0_3_3[j] * IP_0_4_57) + (IP_J_0_3_6[j] * IP_0_4_40) + (IP_J_0_3_8[j] * IP_0_4_71) + (IP_J_0_3_4[j] * IP_0_4_33) + (IP_J_0_3_5[j] * IP_0_4_5) + (IP_J_0_3_0[j] * IP_0_4_12) + (IP_J_0_3_1[j] * IP_0_4_8);
      IP_J_0_6_7[j] = (IP_J_0_3_2[j] * IP_0_4_77) + (IP_J_0_3_7[j] * IP_0_4_38) + (IP_J_0_3_3[j] * IP_0_4_67) + (IP_J_0_3_6[j] * IP_0_4_32) + (IP_J_0_3_8[j] * IP_0_4_43) + (IP_J_0_3_4[j] * IP_0_4_61) + (IP_J_0_3_5[j] * IP_0_4_53) + (IP_J_0_3_0[j] * IP_0_4_17) + (IP_J_0_3_1[j] * IP_0_4_21);
      IP_J_0_6_8[j] = (IP_J_0_3_2[j] * IP_0_4_59) + (IP_J_0_3_7[j] * IP_0_4_44) + (IP_J_0_3_3[j] * IP_0_4_20) + (IP_J_0_3_6[j] * IP_0_4_13) + (IP_J_0_3_8[j] * IP_0_4_3) + (IP_J_0_3_4[j] * IP_0_4_69) + (IP_J_0_3_5[j] * IP_0_4_30) + (IP_J_0_3_0[j] * IP_0_4_47) + (IP_J_0_3_1[j] * IP_0_4_36);
      
    }
    
    for (int j  = 0; j < 30; j += 1)
    {
      
      for (int k  = 0; k < 10; k += 1)
      {
        A[j][k+10] += ((IP_K_0_3_4_EXP_0_8[k+10] * IP_J_0_6_0[j])) + ((IP_K_0_3_1_EXP_0_4[k+10] * IP_J_0_6_5[j])) + ((IP_K_0_3_0_EXP_0_2[k+10] * IP_J_0_6_8[j]));
        A[j][k] += ((IP_K_0_3_7_EXP_0_3[k] * IP_J_0_6_2[j])) + ((IP_K_0_3_5_EXP_0_6[k] * IP_J_0_6_6[j])) + ((IP_K_0_3_6_EXP_0_1[k] * IP_J_0_6_3[j]));
        A[j][k+20] += ((IP_K_0_3_2_EXP_0_5[k+20] * IP_J_0_6_4[j])) + ((IP_K_0_3_3_EXP_0_7[k+20] * IP_J_0_6_7[j])) + ((IP_K_0_3_8_EXP_0_0[k+20] * IP_J_0_6_1[j]));
        
      }
      
    }
    
  }
  
}
            

        
        void wrap_form_cell_integral_0_otherwise(int start, int end,
                      Mat arg0_0_, int *arg0_0_map0_0, int *arg0_0_map1_0, double *arg1_0, int *arg1_0_map0_0, double *arg2_0, int *arg2_0_map0_0, double *arg3_0, double *arg4_0
                      ) {
  Mat arg0_0_0 = arg0_0_;
  double *arg1_0_vec[12];
    double *arg2_0_vec[30];
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
    arg2_0_vec[0] = arg2_0 + (arg2_0_map0_0[i * 10 + 0])* 3;
    arg2_0_vec[1] = arg2_0 + (arg2_0_map0_0[i * 10 + 1])* 3;
    arg2_0_vec[2] = arg2_0 + (arg2_0_map0_0[i * 10 + 2])* 3;
    arg2_0_vec[3] = arg2_0 + (arg2_0_map0_0[i * 10 + 3])* 3;
    arg2_0_vec[4] = arg2_0 + (arg2_0_map0_0[i * 10 + 4])* 3;
    arg2_0_vec[5] = arg2_0 + (arg2_0_map0_0[i * 10 + 5])* 3;
    arg2_0_vec[6] = arg2_0 + (arg2_0_map0_0[i * 10 + 6])* 3;
    arg2_0_vec[7] = arg2_0 + (arg2_0_map0_0[i * 10 + 7])* 3;
    arg2_0_vec[8] = arg2_0 + (arg2_0_map0_0[i * 10 + 8])* 3;
    arg2_0_vec[9] = arg2_0 + (arg2_0_map0_0[i * 10 + 9])* 3;
    arg2_0_vec[10] = arg2_0 + (arg2_0_map0_0[i * 10 + 0])* 3 + 1;
    arg2_0_vec[11] = arg2_0 + (arg2_0_map0_0[i * 10 + 1])* 3 + 1;
    arg2_0_vec[12] = arg2_0 + (arg2_0_map0_0[i * 10 + 2])* 3 + 1;
    arg2_0_vec[13] = arg2_0 + (arg2_0_map0_0[i * 10 + 3])* 3 + 1;
    arg2_0_vec[14] = arg2_0 + (arg2_0_map0_0[i * 10 + 4])* 3 + 1;
    arg2_0_vec[15] = arg2_0 + (arg2_0_map0_0[i * 10 + 5])* 3 + 1;
    arg2_0_vec[16] = arg2_0 + (arg2_0_map0_0[i * 10 + 6])* 3 + 1;
    arg2_0_vec[17] = arg2_0 + (arg2_0_map0_0[i * 10 + 7])* 3 + 1;
    arg2_0_vec[18] = arg2_0 + (arg2_0_map0_0[i * 10 + 8])* 3 + 1;
    arg2_0_vec[19] = arg2_0 + (arg2_0_map0_0[i * 10 + 9])* 3 + 1;
    arg2_0_vec[20] = arg2_0 + (arg2_0_map0_0[i * 10 + 0])* 3 + 2;
    arg2_0_vec[21] = arg2_0 + (arg2_0_map0_0[i * 10 + 1])* 3 + 2;
    arg2_0_vec[22] = arg2_0 + (arg2_0_map0_0[i * 10 + 2])* 3 + 2;
    arg2_0_vec[23] = arg2_0 + (arg2_0_map0_0[i * 10 + 3])* 3 + 2;
    arg2_0_vec[24] = arg2_0 + (arg2_0_map0_0[i * 10 + 4])* 3 + 2;
    arg2_0_vec[25] = arg2_0 + (arg2_0_map0_0[i * 10 + 5])* 3 + 2;
    arg2_0_vec[26] = arg2_0 + (arg2_0_map0_0[i * 10 + 6])* 3 + 2;
    arg2_0_vec[27] = arg2_0 + (arg2_0_map0_0[i * 10 + 7])* 3 + 2;
    arg2_0_vec[28] = arg2_0 + (arg2_0_map0_0[i * 10 + 8])* 3 + 2;
    arg2_0_vec[29] = arg2_0 + (arg2_0_map0_0[i * 10 + 9])* 3 + 2;
    double buffer_arg0_0[30][30]  = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec, arg2_0_vec, arg3_0, arg4_0);
                    double tmp_buffer_arg0_0[30][30] ;
                    for ( int j = 0; j < 10; j++ ) {
                       for ( int k = 0; k < 3; k++ ) {
                          for ( int l = 0; l < 10; l++ ) {
                             for ( int m = 0; m < 3; m++ ) {
                                tmp_buffer_arg0_0[3*j + k][3*l + m] = buffer_arg0_0[j + 10*k][l + 10*m];
                             }
                          }
                       }
                    }
    MatSetValuesBlockedLocal(arg0_0_0, 10, arg0_0_map0_0 + i * 10,
                                             10, arg0_0_map1_0 + i * 10,
                                             (const PetscScalar *)tmp_buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        