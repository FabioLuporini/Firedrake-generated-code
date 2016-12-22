
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


void form_cell_integral_0_otherwise(double A[20][20], double **coordinate_dofs)
{
  // Number of operations (multiply-add pairs) for Jacobian data:      2
  // Number of operations (multiply-add pairs) for geometry tensor:    0
  // Number of operations (multiply-add pairs) for tensor contraction: 140
  // Total number of operations (multiply-add pairs):                  142
  
  // Compute Jacobian
  double J[9];
  compute_jacobian_tetrahedron_3d(J, coordinate_dofs);
  
  // Compute Jacobian inverse and determinant
  double K[9];
  double detJ;
  compute_jacobian_inverse_tetrahedron_3d(K, detJ, J);
  
  const double det = fabs(detJ);
  
  // Compute geometry tensor
  const double G0_ = det;
  
  // Compute element tensor
  A[0][0] = 0.000595238095238099*G0_;
  A[0][1] = 7.4404761904762e-05*G0_;
  A[0][2] = 7.4404761904762e-05*G0_;
  A[0][3] = 7.44047619047623e-05*G0_;
  A[0][4] = 0.000111607142857144*G0_;
  A[0][5] = 0.000111607142857143*G0_;
  A[0][6] = 0.000111607142857143*G0_;
  A[0][7] = 0.000111607142857142*G0_;
  A[0][8] = 0.000111607142857143*G0_;
  A[0][9] = 0.000111607142857143*G0_;
  A[0][10] = -0.000446428571428571*G0_;
  A[0][11] = 0.000223214285714285*G0_;
  A[0][12] = -0.000446428571428571*G0_;
  A[0][13] = 0.000223214285714286*G0_;
  A[0][14] = -0.000446428571428572*G0_;
  A[0][15] = 0.000223214285714285*G0_;
  A[0][16] = 0.00133928571428571*G0_;
  A[0][17] = 0.000669642857142859*G0_;
  A[0][18] = 0.00066964285714286*G0_;
  A[0][19] = 0.000669642857142859*G0_;
  A[1][0] = 7.44047619047621e-05*G0_;
  A[1][1] = 0.000595238095238098*G0_;
  A[1][2] = 7.44047619047615e-05*G0_;
  A[1][3] = 7.44047619047618e-05*G0_;
  A[1][4] = 0.000111607142857143*G0_;
  A[1][5] = 0.000111607142857143*G0_;
  A[1][6] = -0.000446428571428571*G0_;
  A[1][7] = 0.000223214285714286*G0_;
  A[1][8] = -0.000446428571428572*G0_;
  A[1][9] = 0.000223214285714285*G0_;
  A[1][10] = 0.000111607142857142*G0_;
  A[1][11] = 0.000111607142857143*G0_;
  A[1][12] = 0.000111607142857144*G0_;
  A[1][13] = 0.000111607142857143*G0_;
  A[1][14] = 0.000223214285714285*G0_;
  A[1][15] = -0.000446428571428571*G0_;
  A[1][16] = 0.000669642857142856*G0_;
  A[1][17] = 0.00133928571428571*G0_;
  A[1][18] = 0.000669642857142856*G0_;
  A[1][19] = 0.000669642857142855*G0_;
  A[2][0] = 7.44047619047619e-05*G0_;
  A[2][1] = 7.44047619047615e-05*G0_;
  A[2][2] = 0.000595238095238097*G0_;
  A[2][3] = 7.44047619047617e-05*G0_;
  A[2][4] = -0.000446428571428572*G0_;
  A[2][5] = 0.000223214285714287*G0_;
  A[2][6] = 0.000111607142857143*G0_;
  A[2][7] = 0.000111607142857143*G0_;
  A[2][8] = 0.000223214285714286*G0_;
  A[2][9] = -0.000446428571428572*G0_;
  A[2][10] = 0.000111607142857142*G0_;
  A[2][11] = 0.000111607142857143*G0_;
  A[2][12] = 0.000223214285714286*G0_;
  A[2][13] = -0.000446428571428571*G0_;
  A[2][14] = 0.000111607142857143*G0_;
  A[2][15] = 0.000111607142857143*G0_;
  A[2][16] = 0.000669642857142853*G0_;
  A[2][17] = 0.000669642857142853*G0_;
  A[2][18] = 0.00133928571428571*G0_;
  A[2][19] = 0.000669642857142854*G0_;
  A[3][0] = 7.44047619047623e-05*G0_;
  A[3][1] = 7.44047619047618e-05*G0_;
  A[3][2] = 7.44047619047617e-05*G0_;
  A[3][3] = 0.000595238095238098*G0_;
  A[3][4] = 0.000223214285714286*G0_;
  A[3][5] = -0.000446428571428571*G0_;
  A[3][6] = 0.000223214285714286*G0_;
  A[3][7] = -0.000446428571428571*G0_;
  A[3][8] = 0.000111607142857142*G0_;
  A[3][9] = 0.000111607142857143*G0_;
  A[3][10] = 0.000223214285714285*G0_;
  A[3][11] = -0.000446428571428572*G0_;
  A[3][12] = 0.000111607142857143*G0_;
  A[3][13] = 0.000111607142857143*G0_;
  A[3][14] = 0.000111607142857143*G0_;
  A[3][15] = 0.000111607142857143*G0_;
  A[3][16] = 0.000669642857142858*G0_;
  A[3][17] = 0.000669642857142858*G0_;
  A[3][18] = 0.000669642857142858*G0_;
  A[3][19] = 0.00133928571428571*G0_;
  A[4][0] = 0.000111607142857143*G0_;
  A[4][1] = 0.000111607142857143*G0_;
  A[4][2] = -0.000446428571428572*G0_;
  A[4][3] = 0.000223214285714286*G0_;
  A[4][4] = 0.00401785714285714*G0_;
  A[4][5] = -0.00200892857142857*G0_;
  A[4][6] = 0.0;
  A[4][7] = -0.00100446428571429*G0_;
  A[4][8] = -0.00100446428571429*G0_;
  A[4][9] = 0.00200892857142857*G0_;
  A[4][10] = 0.0;
  A[4][11] = -0.00100446428571429*G0_;
  A[4][12] = -0.00100446428571429*G0_;
  A[4][13] = 0.00200892857142857*G0_;
  A[4][14] = 0.0;
  A[4][15] = 0.0;
  A[4][16] = 0.0;
  A[4][17] = 0.0;
  A[4][18] = -0.00200892857142857*G0_;
  A[4][19] = 0.0;
  A[5][0] = 0.000111607142857143*G0_;
  A[5][1] = 0.000111607142857143*G0_;
  A[5][2] = 0.000223214285714287*G0_;
  A[5][3] = -0.000446428571428571*G0_;
  A[5][4] = -0.00200892857142857*G0_;
  A[5][5] = 0.00401785714285715*G0_;
  A[5][6] = -0.00100446428571429*G0_;
  A[5][7] = 0.00200892857142857*G0_;
  A[5][8] = 0.0;
  A[5][9] = -0.00100446428571429*G0_;
  A[5][10] = -0.00100446428571429*G0_;
  A[5][11] = 0.00200892857142857*G0_;
  A[5][12] = 0.0;
  A[5][13] = -0.00100446428571429*G0_;
  A[5][14] = 0.0;
  A[5][15] = 0.0;
  A[5][16] = 0.0;
  A[5][17] = 0.0;
  A[5][18] = 0.0;
  A[5][19] = -0.00200892857142857*G0_;
  A[6][0] = 0.000111607142857143*G0_;
  A[6][1] = -0.000446428571428571*G0_;
  A[6][2] = 0.000111607142857143*G0_;
  A[6][3] = 0.000223214285714286*G0_;
  A[6][4] = 0.0;
  A[6][5] = -0.00100446428571429*G0_;
  A[6][6] = 0.00401785714285714*G0_;
  A[6][7] = -0.00200892857142857*G0_;
  A[6][8] = 0.00200892857142857*G0_;
  A[6][9] = -0.00100446428571429*G0_;
  A[6][10] = 0.0;
  A[6][11] = -0.00100446428571429*G0_;
  A[6][12] = 0.0;
  A[6][13] = 0.0;
  A[6][14] = -0.00100446428571429*G0_;
  A[6][15] = 0.00200892857142857*G0_;
  A[6][16] = 0.0;
  A[6][17] = -0.00200892857142857*G0_;
  A[6][18] = 0.0;
  A[6][19] = 0.0;
  A[7][0] = 0.000111607142857142*G0_;
  A[7][1] = 0.000223214285714286*G0_;
  A[7][2] = 0.000111607142857143*G0_;
  A[7][3] = -0.000446428571428571*G0_;
  A[7][4] = -0.00100446428571429*G0_;
  A[7][5] = 0.00200892857142857*G0_;
  A[7][6] = -0.00200892857142857*G0_;
  A[7][7] = 0.00401785714285715*G0_;
  A[7][8] = -0.00100446428571429*G0_;
  A[7][9] = 0.0;
  A[7][10] = -0.00100446428571429*G0_;
  A[7][11] = 0.00200892857142857*G0_;
  A[7][12] = 0.0;
  A[7][13] = 0.0;
  A[7][14] = 0.0;
  A[7][15] = -0.00100446428571429*G0_;
  A[7][16] = 0.0;
  A[7][17] = 0.0;
  A[7][18] = 0.0;
  A[7][19] = -0.00200892857142858*G0_;
  A[8][0] = 0.000111607142857143*G0_;
  A[8][1] = -0.000446428571428572*G0_;
  A[8][2] = 0.000223214285714286*G0_;
  A[8][3] = 0.000111607142857142*G0_;
  A[8][4] = -0.00100446428571429*G0_;
  A[8][5] = 0.0;
  A[8][6] = 0.00200892857142857*G0_;
  A[8][7] = -0.00100446428571429*G0_;
  A[8][8] = 0.00401785714285714*G0_;
  A[8][9] = -0.00200892857142857*G0_;
  A[8][10] = 0.0;
  A[8][11] = 0.0;
  A[8][12] = 0.0;
  A[8][13] = -0.00100446428571429*G0_;
  A[8][14] = -0.00100446428571429*G0_;
  A[8][15] = 0.00200892857142857*G0_;
  A[8][16] = 0.0;
  A[8][17] = -0.00200892857142857*G0_;
  A[8][18] = 0.0;
  A[8][19] = 0.0;
  A[9][0] = 0.000111607142857143*G0_;
  A[9][1] = 0.000223214285714285*G0_;
  A[9][2] = -0.000446428571428572*G0_;
  A[9][3] = 0.000111607142857143*G0_;
  A[9][4] = 0.00200892857142857*G0_;
  A[9][5] = -0.00100446428571429*G0_;
  A[9][6] = -0.00100446428571429*G0_;
  A[9][7] = 0.0;
  A[9][8] = -0.00200892857142857*G0_;
  A[9][9] = 0.00401785714285714*G0_;
  A[9][10] = 0.0;
  A[9][11] = 0.0;
  A[9][12] = -0.00100446428571429*G0_;
  A[9][13] = 0.00200892857142857*G0_;
  A[9][14] = 0.0;
  A[9][15] = -0.00100446428571429*G0_;
  A[9][16] = 0.0;
  A[9][17] = 0.0;
  A[9][18] = -0.00200892857142857*G0_;
  A[9][19] = 0.0;
  A[10][0] = -0.000446428571428571*G0_;
  A[10][1] = 0.000111607142857142*G0_;
  A[10][2] = 0.000111607142857142*G0_;
  A[10][3] = 0.000223214285714285*G0_;
  A[10][4] = 0.0;
  A[10][5] = -0.00100446428571429*G0_;
  A[10][6] = 0.0;
  A[10][7] = -0.00100446428571429*G0_;
  A[10][8] = 0.0;
  A[10][9] = 0.0;
  A[10][10] = 0.00401785714285714*G0_;
  A[10][11] = -0.00200892857142857*G0_;
  A[10][12] = 0.00200892857142857*G0_;
  A[10][13] = -0.00100446428571429*G0_;
  A[10][14] = 0.00200892857142857*G0_;
  A[10][15] = -0.00100446428571429*G0_;
  A[10][16] = -0.00200892857142858*G0_;
  A[10][17] = 0.0;
  A[10][18] = 0.0;
  A[10][19] = 0.0;
  A[11][0] = 0.000223214285714285*G0_;
  A[11][1] = 0.000111607142857143*G0_;
  A[11][2] = 0.000111607142857143*G0_;
  A[11][3] = -0.000446428571428572*G0_;
  A[11][4] = -0.00100446428571429*G0_;
  A[11][5] = 0.00200892857142857*G0_;
  A[11][6] = -0.00100446428571429*G0_;
  A[11][7] = 0.00200892857142857*G0_;
  A[11][8] = 0.0;
  A[11][9] = 0.0;
  A[11][10] = -0.00200892857142857*G0_;
  A[11][11] = 0.00401785714285715*G0_;
  A[11][12] = -0.00100446428571429*G0_;
  A[11][13] = 0.0;
  A[11][14] = -0.00100446428571429*G0_;
  A[11][15] = 0.0;
  A[11][16] = 0.0;
  A[11][17] = 0.0;
  A[11][18] = 0.0;
  A[11][19] = -0.00200892857142858*G0_;
  A[12][0] = -0.000446428571428571*G0_;
  A[12][1] = 0.000111607142857144*G0_;
  A[12][2] = 0.000223214285714286*G0_;
  A[12][3] = 0.000111607142857143*G0_;
  A[12][4] = -0.00100446428571429*G0_;
  A[12][5] = 0.0;
  A[12][6] = 0.0;
  A[12][7] = 0.0;
  A[12][8] = 0.0;
  A[12][9] = -0.00100446428571429*G0_;
  A[12][10] = 0.00200892857142857*G0_;
  A[12][11] = -0.00100446428571429*G0_;
  A[12][12] = 0.00401785714285715*G0_;
  A[12][13] = -0.00200892857142857*G0_;
  A[12][14] = 0.00200892857142857*G0_;
  A[12][15] = -0.00100446428571429*G0_;
  A[12][16] = -0.00200892857142857*G0_;
  A[12][17] = 0.0;
  A[12][18] = 0.0;
  A[12][19] = 0.0;
  A[13][0] = 0.000223214285714286*G0_;
  A[13][1] = 0.000111607142857143*G0_;
  A[13][2] = -0.000446428571428571*G0_;
  A[13][3] = 0.000111607142857143*G0_;
  A[13][4] = 0.00200892857142857*G0_;
  A[13][5] = -0.00100446428571429*G0_;
  A[13][6] = 0.0;
  A[13][7] = 0.0;
  A[13][8] = -0.00100446428571429*G0_;
  A[13][9] = 0.00200892857142857*G0_;
  A[13][10] = -0.00100446428571429*G0_;
  A[13][11] = 0.0;
  A[13][12] = -0.00200892857142857*G0_;
  A[13][13] = 0.00401785714285715*G0_;
  A[13][14] = -0.00100446428571429*G0_;
  A[13][15] = 0.0;
  A[13][16] = 0.0;
  A[13][17] = 0.0;
  A[13][18] = -0.00200892857142857*G0_;
  A[13][19] = 0.0;
  A[14][0] = -0.000446428571428572*G0_;
  A[14][1] = 0.000223214285714285*G0_;
  A[14][2] = 0.000111607142857143*G0_;
  A[14][3] = 0.000111607142857143*G0_;
  A[14][4] = 0.0;
  A[14][5] = 0.0;
  A[14][6] = -0.00100446428571429*G0_;
  A[14][7] = 0.0;
  A[14][8] = -0.00100446428571429*G0_;
  A[14][9] = 0.0;
  A[14][10] = 0.00200892857142857*G0_;
  A[14][11] = -0.00100446428571429*G0_;
  A[14][12] = 0.00200892857142857*G0_;
  A[14][13] = -0.00100446428571429*G0_;
  A[14][14] = 0.00401785714285715*G0_;
  A[14][15] = -0.00200892857142857*G0_;
  A[14][16] = -0.00200892857142858*G0_;
  A[14][17] = 0.0;
  A[14][18] = 0.0;
  A[14][19] = 0.0;
  A[15][0] = 0.000223214285714285*G0_;
  A[15][1] = -0.000446428571428571*G0_;
  A[15][2] = 0.000111607142857143*G0_;
  A[15][3] = 0.000111607142857143*G0_;
  A[15][4] = 0.0;
  A[15][5] = 0.0;
  A[15][6] = 0.00200892857142857*G0_;
  A[15][7] = -0.00100446428571429*G0_;
  A[15][8] = 0.00200892857142857*G0_;
  A[15][9] = -0.00100446428571429*G0_;
  A[15][10] = -0.00100446428571429*G0_;
  A[15][11] = 0.0;
  A[15][12] = -0.00100446428571429*G0_;
  A[15][13] = 0.0;
  A[15][14] = -0.00200892857142857*G0_;
  A[15][15] = 0.00401785714285715*G0_;
  A[15][16] = 0.0;
  A[15][17] = -0.00200892857142857*G0_;
  A[15][18] = 0.0;
  A[15][19] = 0.0;
  A[16][0] = 0.00133928571428571*G0_;
  A[16][1] = 0.000669642857142856*G0_;
  A[16][2] = 0.000669642857142854*G0_;
  A[16][3] = 0.000669642857142858*G0_;
  A[16][4] = 0.0;
  A[16][5] = 0.0;
  A[16][6] = 0.0;
  A[16][7] = 0.0;
  A[16][8] = 0.0;
  A[16][9] = 0.0;
  A[16][10] = -0.00200892857142858*G0_;
  A[16][11] = 0.0;
  A[16][12] = -0.00200892857142857*G0_;
  A[16][13] = 0.0;
  A[16][14] = -0.00200892857142858*G0_;
  A[16][15] = 0.0;
  A[16][16] = 0.0160714285714286*G0_;
  A[16][17] = 0.00803571428571428*G0_;
  A[16][18] = 0.0080357142857143*G0_;
  A[16][19] = 0.00803571428571429*G0_;
  A[17][0] = 0.000669642857142859*G0_;
  A[17][1] = 0.00133928571428571*G0_;
  A[17][2] = 0.000669642857142853*G0_;
  A[17][3] = 0.000669642857142858*G0_;
  A[17][4] = 0.0;
  A[17][5] = 0.0;
  A[17][6] = -0.00200892857142857*G0_;
  A[17][7] = 0.0;
  A[17][8] = -0.00200892857142857*G0_;
  A[17][9] = 0.0;
  A[17][10] = 0.0;
  A[17][11] = 0.0;
  A[17][12] = 0.0;
  A[17][13] = 0.0;
  A[17][14] = 0.0;
  A[17][15] = -0.00200892857142857*G0_;
  A[17][16] = 0.00803571428571428*G0_;
  A[17][17] = 0.0160714285714286*G0_;
  A[17][18] = 0.00803571428571429*G0_;
  A[17][19] = 0.00803571428571428*G0_;
  A[18][0] = 0.00066964285714286*G0_;
  A[18][1] = 0.000669642857142856*G0_;
  A[18][2] = 0.00133928571428571*G0_;
  A[18][3] = 0.000669642857142858*G0_;
  A[18][4] = -0.00200892857142857*G0_;
  A[18][5] = 0.0;
  A[18][6] = 0.0;
  A[18][7] = 0.0;
  A[18][8] = 0.0;
  A[18][9] = -0.00200892857142857*G0_;
  A[18][10] = 0.0;
  A[18][11] = 0.0;
  A[18][12] = 0.0;
  A[18][13] = -0.00200892857142857*G0_;
  A[18][14] = 0.0;
  A[18][15] = 0.0;
  A[18][16] = 0.0080357142857143*G0_;
  A[18][17] = 0.00803571428571429*G0_;
  A[18][18] = 0.0160714285714286*G0_;
  A[18][19] = 0.00803571428571429*G0_;
  A[19][0] = 0.000669642857142859*G0_;
  A[19][1] = 0.000669642857142855*G0_;
  A[19][2] = 0.000669642857142854*G0_;
  A[19][3] = 0.00133928571428571*G0_;
  A[19][4] = 0.0;
  A[19][5] = -0.00200892857142857*G0_;
  A[19][6] = 0.0;
  A[19][7] = -0.00200892857142857*G0_;
  A[19][8] = 0.0;
  A[19][9] = 0.0;
  A[19][10] = 0.0;
  A[19][11] = -0.00200892857142858*G0_;
  A[19][12] = 0.0;
  A[19][13] = 0.0;
  A[19][14] = 0.0;
  A[19][15] = 0.0;
  A[19][16] = 0.00803571428571429*G0_;
  A[19][17] = 0.00803571428571428*G0_;
  A[19][18] = 0.00803571428571429*G0_;
  A[19][19] = 0.0160714285714286*G0_;
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
    double buffer_arg0_0[20][20] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec);
    MatSetValuesLocal(arg0_0_0, 20, arg0_0_map0_0 + i * 20,
                                             20, arg0_0_map1_0 + i * 20,
                                             (const PetscScalar *)buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        