
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[12][12] , double** coordinate_dofs , double** w0 , double* c1 , double* c2 )
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
  
  
  static const double W1  = 0.166666666666667;
  static const double FE0_C2_D001[1][4]  = {{-1.0, 0.0, 0.0, 1.0}};
  static const double FE0_C1_D010[1][4]  = {{-1.0, 0.0, 1.0, 0.0}};
  static const double FE0_C1_D100[1][4]  = {{-1.0, 1.0, 0.0, 0.0}};
  static const double FE1[1][1]  = {{1.0}};
  double F0  = 0.0;
  double F1  = 0.0;
  double F2  = 0.0;
  double F3  = 0.0;
  double F4  = 0.0;
  double F5  = 0.0;
  double F6  = 0.0;
  double F7  = 0.0;
  double F8  = 0.0;
  double F9  = 0.0;
  double F10  = 0.0;
  
  
  for (int r  = 0; r < 1; r += 1)
  {
    F9 += (w1[0][0] * FE1[0][0]);
    F10 += (w2[0][0] * FE1[0][0]);
    
  }
  
  
  for (int r  = 0; r < 2; r += 1)
  {
    F3 += (w0[r+4][0] * FE0_C1_D100[0][r]);
    F6 += (w0[r+8][0] * FE0_C1_D100[0][r]);
    F0 += (w0[r][0] * FE0_C1_D100[0][r]);
    
  }
  
  
  for (int r  = 0; r < 3; r += 1)
  {
    F4 += (w0[r+4][0] * FE0_C1_D010[0][r]);
    F7 += (w0[r+8][0] * FE0_C1_D010[0][r]);
    F1 += (w0[r][0] * FE0_C1_D010[0][r]);
    
  }
  
  
  for (int r  = 0; r < 4; r += 1)
  {
    F5 += (w0[r+4][0] * FE0_C2_D001[0][r]);
    F8 += (w0[r+8][0] * FE0_C2_D001[0][r]);
    F2 += (w0[r][0] * FE0_C2_D001[0][r]);
    
  }
  
  double CONST_0_1_0  = (2 * 0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0))));
  double CONST_0_1_1  = (2 * 0.5 * ((((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + (((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0))));
  double CONST_0_1_2  = ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0);
  double CONST_0_1_3  = ((F9 * 2 * 0.5 * ((((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0)) + -1.0)) + (F10 * ((0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[8] * F5) + (K[5] * F4) + (K[2] * F3))) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[8] * F2) + (K[5] * F1) + (K[2] * F0))) + (((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + -1.0)) + (0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[7] * F8) + (K[4] * F7) + (K[1] * F6))) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + -1.0)) + (0.5 * ((((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0)) + -1.0)))));
  double CONST_0_1_4  = ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0);
  double CONST_0_1_5  = ((K[6] * F5) + (K[3] * F4) + (K[0] * F3));
  double CONST_0_1_6  = ((F9 * 2 * 0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[7] * F8) + (K[4] * F7) + (K[1] * F6))) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + -1.0)) + (F10 * ((0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[8] * F5) + (K[5] * F4) + (K[2] * F3))) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[8] * F2) + (K[5] * F1) + (K[2] * F0))) + (((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + -1.0)) + (0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[7] * F8) + (K[4] * F7) + (K[1] * F6))) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + -1.0)) + (0.5 * ((((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0)) + -1.0)))));
  double CONST_0_1_7  = ((K[6] * F8) + (K[3] * F7) + (K[0] * F6));
  double CONST_0_1_8  = ((K[8] * F2) + (K[5] * F1) + (K[2] * F0));
  double CONST_0_1_9  = ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0);
  double CONST_0_1_10  = (2 * 0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0))));
  double CONST_0_1_11  = ((K[8] * F5) + (K[5] * F4) + (K[2] * F3));
  double CONST_0_1_12  = ((K[7] * F8) + (K[4] * F7) + (K[1] * F6));
  double CONST_0_1_13  = ((K[7] * F2) + (K[4] * F1) + (K[1] * F0));
  double CONST_0_1_14  = ((F9 * 2 * 0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[8] * F5) + (K[5] * F4) + (K[2] * F3))) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[8] * F2) + (K[5] * F1) + (K[2] * F0))) + (((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + -1.0)) + (F10 * ((0.5 * ((((K[8] * F5) + (K[5] * F4) + (K[2] * F3)) * ((K[8] * F5) + (K[5] * F4) + (K[2] * F3))) + (((K[8] * F2) + (K[5] * F1) + (K[2] * F0)) * ((K[8] * F2) + (K[5] * F1) + (K[2] * F0))) + (((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0) * ((K[8] * F8) + (K[5] * F7) + (K[2] * F6) + 1.0)) + -1.0)) + (0.5 * ((((K[7] * F8) + (K[4] * F7) + (K[1] * F6)) * ((K[7] * F8) + (K[4] * F7) + (K[1] * F6))) + (((K[7] * F2) + (K[4] * F1) + (K[1] * F0)) * ((K[7] * F2) + (K[4] * F1) + (K[1] * F0))) + (((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0) * ((K[7] * F5) + (K[4] * F4) + (K[1] * F3) + 1.0)) + -1.0)) + (0.5 * ((((K[6] * F8) + (K[3] * F7) + (K[0] * F6)) * ((K[6] * F8) + (K[3] * F7) + (K[0] * F6))) + (((K[6] * F5) + (K[3] * F4) + (K[0] * F3)) * ((K[6] * F5) + (K[3] * F4) + (K[0] * F3))) + (((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0) * ((K[6] * F2) + (K[3] * F1) + (K[0] * F0) + 1.0)) + -1.0)))));
  double J_0_3_0[12]  = {0.0};
  double J_0_3_1[12]  = {0.0};
  double J_0_3_2[12]  = {0.0};
  double J_0_3_3[12]  = {0.0};
  double J_0_3_4[12]  = {0.0};
  double J_0_3_5[12]  = {0.0};
  double J_0_3_6[12]  = {0.0};
  double J_0_3_7[12]  = {0.0};
  double J_0_3_8[12]  = {0.0};
  double K_0_3_0[12]  = {0.0};
  double K_0_3_1[12]  = {0.0};
  double K_0_3_2[12]  = {0.0};
  double K_0_3_3[12]  = {0.0};
  double K_0_3_4[12]  = {0.0};
  double K_0_3_5[12]  = {0.0};
  double K_0_3_6[12]  = {0.0};
  double K_0_3_7[12]  = {0.0};
  double K_0_3_8[12]  = {0.0};
  double CONST_0_4_0  = ((CONST_0_1_11 * CONST_0_1_7 * F9) + (CONST_0_1_5 * F10 * CONST_0_1_4));
  double CONST_0_4_1  = ((CONST_0_1_11 * CONST_0_1_12 * F9) + (CONST_0_1_9 * F10 * CONST_0_1_4));
  double CONST_0_4_2  = ((CONST_0_1_7 * CONST_0_1_9 * F9) + (CONST_0_1_12 * F10 * CONST_0_1_5));
  double CONST_0_4_3  = ((CONST_0_1_2 * CONST_0_1_12 * F9) + (CONST_0_1_13 * F10 * CONST_0_1_7));
  double CONST_0_4_4  = ((CONST_0_1_5 * CONST_0_1_8 * F9) + (CONST_0_1_11 * F10 * CONST_0_1_2));
  double CONST_0_4_5  = ((CONST_0_1_13 * CONST_0_1_5 * F9) + (CONST_0_1_2 * F10 * CONST_0_1_9));
  double CONST_0_4_6  = ((CONST_0_1_2 * CONST_0_1_11 * F9) + (CONST_0_1_8 * F10 * CONST_0_1_5));
  double CONST_0_4_7  = ((CONST_0_1_13 * CONST_0_1_11 * F9) + (CONST_0_1_8 * F10 * CONST_0_1_9));
  double CONST_0_4_8  = ((CONST_0_1_12 * CONST_0_1_5 * F9) + (CONST_0_1_7 * F10 * CONST_0_1_9));
  double CONST_0_4_9  = ((CONST_0_1_4 * CONST_0_1_2 * F9) + (CONST_0_1_7 * F10 * CONST_0_1_8));
  double CONST_0_4_10  = ((CONST_0_1_11 * CONST_0_1_5 * F9) + (CONST_0_1_10 * F9) + (CONST_0_1_5 * CONST_0_1_11 * F10));
  double CONST_0_4_11  = ((F9 * ((CONST_0_1_4 * CONST_0_1_11) + (CONST_0_1_12 * CONST_0_1_9) + (2.0 * CONST_0_1_7 * CONST_0_1_5))) + (CONST_0_1_7 * F10 * CONST_0_1_5));
  double CONST_0_4_12  = ((CONST_0_1_5 * CONST_0_1_13 * F9) + (CONST_0_1_9 * F10 * CONST_0_1_2));
  double CONST_0_4_13  = ((CONST_0_1_8 * CONST_0_1_8 * F9) + (CONST_0_1_13 * CONST_0_1_13 * F9) + CONST_0_1_3 + (CONST_0_1_2 * CONST_0_1_2 * ((2.0 * F9) + F10)));
  double CONST_0_4_14  = ((F9 * ((CONST_0_1_9 * CONST_0_1_12) + (CONST_0_1_5 * CONST_0_1_7) + (2.0 * CONST_0_1_11 * CONST_0_1_4))) + (CONST_0_1_11 * F10 * CONST_0_1_4));
  double CONST_0_4_15  = ((CONST_0_1_13 * CONST_0_1_4 * F9) + (CONST_0_1_8 * F10 * CONST_0_1_12));
  double CONST_0_4_16  = ((CONST_0_1_13 * CONST_0_1_13 * F9) + (CONST_0_1_2 * CONST_0_1_2 * F9) + (CONST_0_1_8 * CONST_0_1_8 * ((2.0 * F9) + F10)) + CONST_0_1_14);
  double CONST_0_4_17  = ((CONST_0_1_8 * CONST_0_1_12 * F9) + (CONST_0_1_13 * F10 * CONST_0_1_4));
  double CONST_0_4_18  = ((F9 * ((CONST_0_1_13 * CONST_0_1_12) + (CONST_0_1_2 * CONST_0_1_7) + (2.0 * CONST_0_1_8 * CONST_0_1_4))) + (CONST_0_1_8 * F10 * CONST_0_1_4));
  double CONST_0_4_19  = ((F9 * ((CONST_0_1_4 * CONST_0_1_8) + (CONST_0_1_12 * CONST_0_1_13) + (2.0 * CONST_0_1_7 * CONST_0_1_2))) + (CONST_0_1_7 * F10 * CONST_0_1_2));
  double CONST_0_4_20  = ((CONST_0_1_13 * CONST_0_1_8 * F9) + (CONST_0_1_1 * F9) + (CONST_0_1_8 * CONST_0_1_13 * F10));
  double CONST_0_4_21  = ((CONST_0_1_9 * CONST_0_1_4 * F9) + (CONST_0_1_11 * F10 * CONST_0_1_12));
  double CONST_0_4_22  = ((CONST_0_1_9 * CONST_0_1_7 * F9) + (CONST_0_1_5 * F10 * CONST_0_1_12));
  double CONST_0_4_23  = ((F9 * ((CONST_0_1_8 * CONST_0_1_11) + (CONST_0_1_13 * CONST_0_1_9) + (2.0 * CONST_0_1_2 * CONST_0_1_5))) + (CONST_0_1_2 * F10 * CONST_0_1_5));
  double CONST_0_4_24  = ((CONST_0_1_13 * CONST_0_1_7 * F9) + (CONST_0_1_2 * F10 * CONST_0_1_12));
  double CONST_0_4_25  = ((F9 * (CONST_0_1_10 + (CONST_0_1_7 * CONST_0_1_4))) + (CONST_0_1_4 * F10 * CONST_0_1_7));
  double CONST_0_4_26  = ((CONST_0_1_7 * CONST_0_1_12 * F9) + (CONST_0_1_0 * F9) + (CONST_0_1_12 * CONST_0_1_7 * F10));
  double CONST_0_4_27  = ((F9 * (CONST_0_1_1 + (CONST_0_1_11 * CONST_0_1_9))) + (CONST_0_1_9 * F10 * CONST_0_1_11));
  double CONST_0_4_28  = ((CONST_0_1_2 * CONST_0_1_4 * F9) + (CONST_0_1_8 * F10 * CONST_0_1_7));
  double CONST_0_4_29  = ((CONST_0_1_2 * CONST_0_1_9 * F9) + (CONST_0_1_13 * F10 * CONST_0_1_5));
  double CONST_0_4_30  = ((CONST_0_1_4 * CONST_0_1_5 * F9) + (CONST_0_1_7 * F10 * CONST_0_1_11));
  double CONST_0_4_31  = ((CONST_0_1_4 * CONST_0_1_7 * F9) + (CONST_0_1_10 * F9) + (CONST_0_1_7 * CONST_0_1_4 * F10));
  double CONST_0_4_32  = ((F9 * ((CONST_0_1_4 * CONST_0_1_11) + (CONST_0_1_7 * CONST_0_1_5) + (2.0 * CONST_0_1_12 * CONST_0_1_9))) + (CONST_0_1_12 * F10 * CONST_0_1_9));
  double CONST_0_4_33  = ((CONST_0_1_7 * CONST_0_1_13 * F9) + (CONST_0_1_12 * F10 * CONST_0_1_2));
  double CONST_0_4_34  = ((CONST_0_1_11 * CONST_0_1_11 * F9) + (CONST_0_1_5 * CONST_0_1_5 * F9) + CONST_0_1_6 + (CONST_0_1_9 * CONST_0_1_9 * ((2.0 * F9) + F10)));
  double CONST_0_4_35  = ((F9 * ((CONST_0_1_9 * CONST_0_1_12) + (CONST_0_1_11 * CONST_0_1_4) + (2.0 * CONST_0_1_5 * CONST_0_1_7))) + (CONST_0_1_5 * F10 * CONST_0_1_7));
  double CONST_0_4_36  = ((CONST_0_1_2 * CONST_0_1_13 * F9) + (CONST_0_1_0 * F9) + (CONST_0_1_13 * CONST_0_1_2 * F10));
  double CONST_0_4_37  = ((CONST_0_1_7 * CONST_0_1_7 * F9) + (CONST_0_1_4 * CONST_0_1_4 * F9) + (CONST_0_1_12 * CONST_0_1_12 * ((2.0 * F9) + F10)) + CONST_0_1_6);
  double CONST_0_4_38  = ((F9 * ((CONST_0_1_9 * CONST_0_1_13) + (CONST_0_1_5 * CONST_0_1_2) + (2.0 * CONST_0_1_11 * CONST_0_1_8))) + (CONST_0_1_11 * F10 * CONST_0_1_8));
  double CONST_0_4_39  = ((CONST_0_1_7 * CONST_0_1_11 * F9) + (CONST_0_1_4 * F10 * CONST_0_1_5));
  double CONST_0_4_40  = ((CONST_0_1_8 * CONST_0_1_5 * F9) + (CONST_0_1_2 * F10 * CONST_0_1_11));
  double CONST_0_4_41  = ((CONST_0_1_4 * CONST_0_1_13 * F9) + (CONST_0_1_12 * F10 * CONST_0_1_8));
  double CONST_0_4_42  = ((F9 * ((CONST_0_1_8 * CONST_0_1_11) + (CONST_0_1_2 * CONST_0_1_5) + (2.0 * CONST_0_1_13 * CONST_0_1_9))) + (CONST_0_1_13 * F10 * CONST_0_1_9));
  double CONST_0_4_43  = ((F9 * ((CONST_0_1_2 * CONST_0_1_7) + (CONST_0_1_8 * CONST_0_1_4) + (2.0 * CONST_0_1_13 * CONST_0_1_12))) + (CONST_0_1_13 * F10 * CONST_0_1_12));
  double CONST_0_4_44  = ((F9 * (CONST_0_1_10 + (CONST_0_1_8 * CONST_0_1_2))) + (CONST_0_1_2 * F10 * CONST_0_1_8));
  double CONST_0_4_45  = ((CONST_0_1_8 * CONST_0_1_8 * F9) + (CONST_0_1_2 * CONST_0_1_2 * F9) + (CONST_0_1_13 * CONST_0_1_13 * ((2.0 * F9) + F10)) + CONST_0_1_6);
  double CONST_0_4_46  = ((CONST_0_1_11 * CONST_0_1_2 * F9) + (CONST_0_1_5 * F10 * CONST_0_1_8));
  double CONST_0_4_47  = ((F9 * ((CONST_0_1_11 * CONST_0_1_8) + (CONST_0_1_9 * CONST_0_1_13) + (2.0 * CONST_0_1_5 * CONST_0_1_2))) + (CONST_0_1_5 * F10 * CONST_0_1_2));
  double CONST_0_4_48  = ((F9 * (CONST_0_1_0 + (CONST_0_1_13 * CONST_0_1_2))) + (CONST_0_1_2 * F10 * CONST_0_1_13));
  double CONST_0_4_49  = ((CONST_0_1_11 * CONST_0_1_13 * F9) + (CONST_0_1_9 * F10 * CONST_0_1_8));
  double CONST_0_4_50  = ((F9 * (CONST_0_1_1 + (CONST_0_1_12 * CONST_0_1_4))) + (CONST_0_1_4 * F10 * CONST_0_1_12));
  double CONST_0_4_51  = ((CONST_0_1_4 * CONST_0_1_12 * F9) + (CONST_0_1_1 * F9) + (CONST_0_1_12 * CONST_0_1_4 * F10));
  double CONST_0_4_52  = ((CONST_0_1_2 * CONST_0_1_8 * F9) + (CONST_0_1_10 * F9) + (CONST_0_1_8 * CONST_0_1_2 * F10));
  double CONST_0_4_53  = ((CONST_0_1_8 * CONST_0_1_13 * F9) + (CONST_0_1_1 * F9) + (CONST_0_1_13 * CONST_0_1_8 * F10));
  double CONST_0_4_54  = ((F9 * ((CONST_0_1_12 * CONST_0_1_13) + (CONST_0_1_7 * CONST_0_1_2) + (2.0 * CONST_0_1_4 * CONST_0_1_8))) + (CONST_0_1_4 * F10 * CONST_0_1_8));
  double CONST_0_4_55  = ((F9 * ((CONST_0_1_4 * CONST_0_1_8) + (CONST_0_1_7 * CONST_0_1_2) + (2.0 * CONST_0_1_12 * CONST_0_1_13))) + (CONST_0_1_12 * F10 * CONST_0_1_13));
  double CONST_0_4_56  = ((CONST_0_1_12 * CONST_0_1_2 * F9) + (CONST_0_1_7 * F10 * CONST_0_1_13));
  double CONST_0_4_57  = ((CONST_0_1_5 * CONST_0_1_4 * F9) + (CONST_0_1_11 * F10 * CONST_0_1_7));
  double CONST_0_4_58  = ((CONST_0_1_12 * CONST_0_1_12 * F9) + (CONST_0_1_4 * CONST_0_1_4 * F9) + (CONST_0_1_7 * CONST_0_1_7 * ((2.0 * F9) + F10)) + CONST_0_1_3);
  double CONST_0_4_59  = ((CONST_0_1_8 * CONST_0_1_7 * F9) + (CONST_0_1_2 * F10 * CONST_0_1_4));
  double CONST_0_4_60  = ((CONST_0_1_5 * CONST_0_1_12 * F9) + (CONST_0_1_9 * F10 * CONST_0_1_7));
  double CONST_0_4_61  = ((CONST_0_1_7 * CONST_0_1_8 * F9) + (CONST_0_1_4 * F10 * CONST_0_1_2));
  double CONST_0_4_62  = ((F9 * ((CONST_0_1_13 * CONST_0_1_12) + (CONST_0_1_8 * CONST_0_1_4) + (2.0 * CONST_0_1_2 * CONST_0_1_7))) + (CONST_0_1_2 * F10 * CONST_0_1_7));
  double CONST_0_4_63  = ((CONST_0_1_9 * CONST_0_1_5 * F9) + (CONST_0_1_0 * F9) + (CONST_0_1_5 * CONST_0_1_9 * F10));
  double CONST_0_4_64  = ((CONST_0_1_12 * CONST_0_1_7 * F9) + (CONST_0_1_0 * F9) + (CONST_0_1_7 * CONST_0_1_12 * F10));
  double CONST_0_4_65  = ((F9 * ((CONST_0_1_11 * CONST_0_1_8) + (CONST_0_1_5 * CONST_0_1_2) + (2.0 * CONST_0_1_9 * CONST_0_1_13))) + (CONST_0_1_9 * F10 * CONST_0_1_13));
  double CONST_0_4_66  = ((CONST_0_1_8 * CONST_0_1_9 * F9) + (CONST_0_1_13 * F10 * CONST_0_1_11));
  double CONST_0_4_67  = ((F9 * ((CONST_0_1_7 * CONST_0_1_5) + (CONST_0_1_12 * CONST_0_1_9) + (2.0 * CONST_0_1_4 * CONST_0_1_11))) + (CONST_0_1_4 * F10 * CONST_0_1_11));
  double CONST_0_4_68  = ((CONST_0_1_12 * CONST_0_1_8 * F9) + (CONST_0_1_4 * F10 * CONST_0_1_13));
  double CONST_0_4_69  = ((F9 * (CONST_0_1_0 + (CONST_0_1_5 * CONST_0_1_9))) + (CONST_0_1_9 * F10 * CONST_0_1_5));
  double CONST_0_4_70  = ((CONST_0_1_5 * CONST_0_1_11 * F9) + (CONST_0_1_10 * F9) + (CONST_0_1_11 * CONST_0_1_5 * F10));
  double CONST_0_4_71  = ((CONST_0_1_9 * CONST_0_1_2 * F9) + (CONST_0_1_5 * F10 * CONST_0_1_13));
  double CONST_0_4_72  = ((CONST_0_1_12 * CONST_0_1_11 * F9) + (CONST_0_1_4 * F10 * CONST_0_1_9));
  double CONST_0_4_73  = ((F9 * ((CONST_0_1_5 * CONST_0_1_7) + (CONST_0_1_11 * CONST_0_1_4) + (2.0 * CONST_0_1_9 * CONST_0_1_12))) + (CONST_0_1_9 * F10 * CONST_0_1_12));
  double CONST_0_4_74  = ((CONST_0_1_9 * CONST_0_1_11 * F9) + (CONST_0_1_1 * F9) + (CONST_0_1_11 * CONST_0_1_9 * F10));
  double CONST_0_4_75  = ((CONST_0_1_9 * CONST_0_1_8 * F9) + (CONST_0_1_11 * F10 * CONST_0_1_13));
  double CONST_0_4_76  = ((CONST_0_1_12 * CONST_0_1_12 * F9) + (CONST_0_1_7 * CONST_0_1_7 * F9) + CONST_0_1_14 + (CONST_0_1_4 * CONST_0_1_4 * ((2.0 * F9) + F10)));
  double CONST_0_4_77  = ((CONST_0_1_11 * CONST_0_1_11 * F9) + (CONST_0_1_9 * CONST_0_1_9 * F9) + (CONST_0_1_5 * CONST_0_1_5 * ((2.0 * F9) + F10)) + CONST_0_1_3);
  double CONST_0_4_78  = ((CONST_0_1_5 * CONST_0_1_5 * F9) + (CONST_0_1_9 * CONST_0_1_9 * F9) + (CONST_0_1_11 * CONST_0_1_11 * ((2.0 * F9) + F10)) + CONST_0_1_14);
  double CONST_0_4_79  = ((CONST_0_1_4 * CONST_0_1_9 * F9) + (CONST_0_1_12 * F10 * CONST_0_1_11));
  double CONST_0_4_80  = ((F9 * ((CONST_0_1_2 * CONST_0_1_5) + (CONST_0_1_13 * CONST_0_1_9) + (2.0 * CONST_0_1_8 * CONST_0_1_11))) + (CONST_0_1_8 * F10 * CONST_0_1_11));
  double K_0_6_0[12]  = {0.0};
  double K_0_6_1[12]  = {0.0};
  double K_0_6_2[12]  = {0.0};
  double K_0_6_3[12]  = {0.0};
  double K_0_6_4[12]  = {0.0};
  double K_0_6_5[12]  = {0.0};
  double K_0_6_6[12]  = {0.0};
  double K_0_6_7[12]  = {0.0};
  double K_0_6_8[12]  = {0.0};
  
  
  for (int k  = 0; k < 4; k += 1)
  {
    J_0_3_0[k] = (K[2] * FE0_C1_D100[0][k]) + (K[5] * FE0_C1_D010[0][k]) + (K[8] * FE0_C2_D001[0][k]);
    J_0_3_1[k+4] = (K[1] * FE0_C1_D100[0][k]) + (K[4] * FE0_C1_D010[0][k]) + (K[7] * FE0_C2_D001[0][k]);
    J_0_3_2[k+8] = (K[2] * FE0_C1_D100[0][k]) + (K[5] * FE0_C1_D010[0][k]) + (K[8] * FE0_C2_D001[0][k]);
    J_0_3_3[k] = (K[0] * FE0_C1_D100[0][k]) + (K[3] * FE0_C1_D010[0][k]) + (K[6] * FE0_C2_D001[0][k]);
    J_0_3_4[k+8] = (K[1] * FE0_C1_D100[0][k]) + (K[4] * FE0_C1_D010[0][k]) + (K[7] * FE0_C2_D001[0][k]);
    J_0_3_5[k+4] = (K[2] * FE0_C1_D100[0][k]) + (K[5] * FE0_C1_D010[0][k]) + (K[8] * FE0_C2_D001[0][k]);
    J_0_3_6[k] = (K[1] * FE0_C1_D100[0][k]) + (K[4] * FE0_C1_D010[0][k]) + (K[7] * FE0_C2_D001[0][k]);
    J_0_3_7[k+8] = (K[0] * FE0_C1_D100[0][k]) + (K[3] * FE0_C1_D010[0][k]) + (K[6] * FE0_C2_D001[0][k]);
    J_0_3_8[k+4] = (K[0] * FE0_C1_D100[0][k]) + (K[3] * FE0_C1_D010[0][k]) + (K[6] * FE0_C2_D001[0][k]);
    K_0_3_0[k] = (FE0_C1_D100[0][k] * K[1]) + (FE0_C1_D010[0][k] * K[4]) + (FE0_C2_D001[0][k] * K[7]);
    K_0_3_1[k] = (FE0_C1_D100[0][k] * K[2]) + (FE0_C1_D010[0][k] * K[5]) + (FE0_C2_D001[0][k] * K[8]);
    K_0_3_2[k+4] = (FE0_C1_D100[0][k] * K[1]) + (FE0_C1_D010[0][k] * K[4]) + (FE0_C2_D001[0][k] * K[7]);
    K_0_3_3[k+8] = (FE0_C1_D100[0][k] * K[2]) + (FE0_C1_D010[0][k] * K[5]) + (FE0_C2_D001[0][k] * K[8]);
    K_0_3_4[k+4] = (FE0_C1_D100[0][k] * K[2]) + (FE0_C1_D010[0][k] * K[5]) + (FE0_C2_D001[0][k] * K[8]);
    K_0_3_5[k] = (FE0_C1_D100[0][k] * K[0]) + (FE0_C1_D010[0][k] * K[3]) + (FE0_C2_D001[0][k] * K[6]);
    K_0_3_6[k+8] = (FE0_C1_D100[0][k] * K[0]) + (FE0_C1_D010[0][k] * K[3]) + (FE0_C2_D001[0][k] * K[6]);
    K_0_3_7[k+8] = (FE0_C1_D100[0][k] * K[1]) + (FE0_C1_D010[0][k] * K[4]) + (FE0_C2_D001[0][k] * K[7]);
    K_0_3_8[k+4] = (FE0_C1_D100[0][k] * K[0]) + (FE0_C1_D010[0][k] * K[3]) + (FE0_C2_D001[0][k] * K[6]);
    
  }
  
  
  for (int k  = 0; k < 12; k += 1)
  {
    K_0_6_0[k] = ((K_0_3_3[k] * CONST_0_4_61) + (K_0_3_0[k] * CONST_0_4_36) + (K_0_3_7[k] * CONST_0_4_33) + (K_0_3_8[k] * CONST_0_4_47) + (K_0_3_6[k] * CONST_0_4_19) + (K_0_3_5[k] * CONST_0_4_13) + (K_0_3_4[k] * CONST_0_4_4) + (K_0_3_2[k] * CONST_0_4_12) + (K_0_3_1[k] * CONST_0_4_52)) * det * W1;
    K_0_6_1[k] = ((K_0_3_3[k] * CONST_0_4_76) + (K_0_3_0[k] * CONST_0_4_17) + (K_0_3_7[k] * CONST_0_4_51) + (K_0_3_8[k] * CONST_0_4_0) + (K_0_3_6[k] * CONST_0_4_31) + (K_0_3_5[k] * CONST_0_4_59) + (K_0_3_4[k] * CONST_0_4_14) + (K_0_3_2[k] * CONST_0_4_1) + (K_0_3_1[k] * CONST_0_4_18)) * det * W1;
    K_0_6_2[k] = ((K_0_3_3[k] * CONST_0_4_50) + (K_0_3_0[k] * CONST_0_4_43) + (K_0_3_7[k] * CONST_0_4_37) + (K_0_3_8[k] * CONST_0_4_22) + (K_0_3_6[k] * CONST_0_4_64) + (K_0_3_5[k] * CONST_0_4_24) + (K_0_3_4[k] * CONST_0_4_21) + (K_0_3_2[k] * CONST_0_4_73) + (K_0_3_1[k] * CONST_0_4_15)) * det * W1;
    K_0_6_3[k] = ((K_0_3_3[k] * CONST_0_4_72) + (K_0_3_0[k] * CONST_0_4_42) + (K_0_3_7[k] * CONST_0_4_32) + (K_0_3_8[k] * CONST_0_4_63) + (K_0_3_6[k] * CONST_0_4_8) + (K_0_3_5[k] * CONST_0_4_5) + (K_0_3_4[k] * CONST_0_4_74) + (K_0_3_2[k] * CONST_0_4_34) + (K_0_3_1[k] * CONST_0_4_7)) * det * W1;
    K_0_6_4[k] = ((K_0_3_3[k] * CONST_0_4_39) + (K_0_3_0[k] * CONST_0_4_29) + (K_0_3_7[k] * CONST_0_4_2) + (K_0_3_8[k] * CONST_0_4_77) + (K_0_3_6[k] * CONST_0_4_11) + (K_0_3_5[k] * CONST_0_4_23) + (K_0_3_4[k] * CONST_0_4_70) + (K_0_3_2[k] * CONST_0_4_69) + (K_0_3_1[k] * CONST_0_4_6)) * det * W1;
    K_0_6_5[k] = ((K_0_3_6[k] * CONST_0_4_58) + (K_0_3_0[k] * CONST_0_4_3) + (K_0_3_7[k] * CONST_0_4_26) + (K_0_3_8[k] * CONST_0_4_35) + (K_0_3_5[k] * CONST_0_4_62) + (K_0_3_3[k] * CONST_0_4_25) + (K_0_3_4[k] * CONST_0_4_57) + (K_0_3_2[k] * CONST_0_4_60) + (K_0_3_1[k] * CONST_0_4_28)) * det * W1;
    K_0_6_6[k] = ((K_0_3_2[k] * CONST_0_4_27) + (K_0_3_0[k] * CONST_0_4_66) + (K_0_3_7[k] * CONST_0_4_79) + (K_0_3_8[k] * CONST_0_4_10) + (K_0_3_6[k] * CONST_0_4_30) + (K_0_3_5[k] * CONST_0_4_40) + (K_0_3_3[k] * CONST_0_4_67) + (K_0_3_4[k] * CONST_0_4_78) + (K_0_3_1[k] * CONST_0_4_80)) * det * W1;
    K_0_6_7[k] = ((K_0_3_3[k] * CONST_0_4_68) + (K_0_3_0[k] * CONST_0_4_45) + (K_0_3_7[k] * CONST_0_4_55) + (K_0_3_8[k] * CONST_0_4_71) + (K_0_3_6[k] * CONST_0_4_56) + (K_0_3_5[k] * CONST_0_4_48) + (K_0_3_4[k] * CONST_0_4_75) + (K_0_3_2[k] * CONST_0_4_65) + (K_0_3_1[k] * CONST_0_4_20)) * det * W1;
    K_0_6_8[k] = ((K_0_3_3[k] * CONST_0_4_54) + (K_0_3_0[k] * CONST_0_4_53) + (K_0_3_7[k] * CONST_0_4_41) + (K_0_3_8[k] * CONST_0_4_46) + (K_0_3_6[k] * CONST_0_4_9) + (K_0_3_5[k] * CONST_0_4_44) + (K_0_3_4[k] * CONST_0_4_38) + (K_0_3_2[k] * CONST_0_4_49) + (K_0_3_1[k] * CONST_0_4_16)) * det * W1;
    
  }
  
  
  
  for (int j  = 0; j < 4; j += 1)
  {
    
    for (int k  = 0; k < 12; k += 1)
    {
      A[j+4][k] += ((K_0_6_4[k] * J_0_3_8[j+4])) + ((K_0_6_3[k] * J_0_3_1[j+4])) + ((K_0_6_6[k] * J_0_3_5[j+4]));
      A[j][k] += ((K_0_6_0[k] * J_0_3_3[j])) + ((K_0_6_7[k] * J_0_3_6[j])) + ((K_0_6_8[k] * J_0_3_0[j]));
      A[j+8][k] += ((K_0_6_1[k] * J_0_3_2[j+8])) + ((K_0_6_5[k] * J_0_3_7[j+8])) + ((K_0_6_2[k] * J_0_3_4[j+8]));
      
    }
    
  }
  
}
            

        
        void wrap_form_cell_integral_0_otherwise(int start, int end,
                      Mat arg0_0_, int *arg0_0_map0_0, int *arg0_0_map1_0, double *arg1_0, int *arg1_0_map0_0, double *arg2_0, int *arg2_0_map0_0, double *arg3_0, double *arg4_0
                      ) {
  Mat arg0_0_0 = arg0_0_;
  double *arg1_0_vec[12];
    double *arg2_0_vec[12];
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
    arg2_0_vec[0] = arg2_0 + (arg2_0_map0_0[i * 4 + 0])* 3;
    arg2_0_vec[1] = arg2_0 + (arg2_0_map0_0[i * 4 + 1])* 3;
    arg2_0_vec[2] = arg2_0 + (arg2_0_map0_0[i * 4 + 2])* 3;
    arg2_0_vec[3] = arg2_0 + (arg2_0_map0_0[i * 4 + 3])* 3;
    arg2_0_vec[4] = arg2_0 + (arg2_0_map0_0[i * 4 + 0])* 3 + 1;
    arg2_0_vec[5] = arg2_0 + (arg2_0_map0_0[i * 4 + 1])* 3 + 1;
    arg2_0_vec[6] = arg2_0 + (arg2_0_map0_0[i * 4 + 2])* 3 + 1;
    arg2_0_vec[7] = arg2_0 + (arg2_0_map0_0[i * 4 + 3])* 3 + 1;
    arg2_0_vec[8] = arg2_0 + (arg2_0_map0_0[i * 4 + 0])* 3 + 2;
    arg2_0_vec[9] = arg2_0 + (arg2_0_map0_0[i * 4 + 1])* 3 + 2;
    arg2_0_vec[10] = arg2_0 + (arg2_0_map0_0[i * 4 + 2])* 3 + 2;
    arg2_0_vec[11] = arg2_0 + (arg2_0_map0_0[i * 4 + 3])* 3 + 2;
    double buffer_arg0_0[12][12] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec, arg2_0_vec, arg3_0, arg4_0);
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
        
        