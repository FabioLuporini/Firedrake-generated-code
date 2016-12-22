
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[30][30] , double** coordinate_dofs )
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
  
  
  static const double W4[4]  = {0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 0.0416666666666667};
  static const double FE1_C1_D100[4][7]  = {{0.447213595499965, 1.34164078649988, 0.552786404500044, 0.552786404500044, -0.552786404500044, -0.552786404500045, -1.78885438199984}, 
  {0.447213595499965, -0.447213595499958, 0.552786404500044, 2.34164078649988, -0.552786404500044, -2.34164078649988, 0.0}, 
  {0.447213595499965, -0.447213595499956, 2.34164078649988, 0.552786404500044, -2.34164078649988, -0.552786404500044, 0.0}, 
  {-1.34164078649987, -0.447213595499956, 0.552786404500044, 0.552786404500045, -0.552786404500044, -0.552786404500045, 1.78885438199982}};
  static const double FE1_C1_D010[4][7]  = {{0.447213595499965, -0.447213595499959, 0.552786404500045, 2.34164078649988, -0.552786404500044, 0.0, -2.34164078649988}, 
  {0.447213595499965, 1.34164078649987, 0.552786404500045, 0.552786404500051, -0.552786404500045, -1.78885438199984, -0.552786404500045}, 
  {0.447213595499965, -0.447213595499958, 2.34164078649988, 0.552786404500048, -2.34164078649988, 0.0, -0.552786404500044}, 
  {-1.34164078649987, -0.447213595499959, 0.552786404500045, 0.55278640450005, -0.552786404500044, 1.78885438199982, -0.552786404500045}};
  static const double FE1_C1_D001[4][7]  = {{0.447213595499965, -0.44721359549996, 0.552786404500046, 2.34164078649988, 0.0, -0.552786404500044, -2.34164078649988}, 
  {0.447213595499964, -0.44721359549996, 2.34164078649988, 0.552786404500049, 0.0, -2.34164078649988, -0.552786404500044}, 
  {0.447213595499965, 1.34164078649987, 0.552786404500042, 0.552786404500052, -1.78885438199984, -0.552786404500043, -0.552786404500045}, 
  {-1.34164078649987, -0.44721359549996, 0.552786404500046, 0.552786404500049, 1.78885438199982, -0.552786404500044, -0.552786404500044}};
  
  const unsigned int nzc18[7] = {20, 23, 24, 25, 27, 28, 29};
  
  const unsigned int nzc21[7] = {0, 3, 4, 5, 7, 8, 9};
  
  const unsigned int nzc14[7] = {10, 13, 14, 15, 17, 18, 19};
  
  const unsigned int nzc22[7] = {0, 2, 4, 6, 7, 8, 9};
  
  const unsigned int nzc19[7] = {20, 22, 24, 26, 27, 28, 29};
  
  const unsigned int nzc15[7] = {10, 12, 14, 16, 17, 18, 19};
  
  const unsigned int nzc23[7] = {0, 1, 5, 6, 7, 8, 9};
  
  const unsigned int nzc16[7] = {10, 11, 15, 16, 17, 18, 19};
  
  const unsigned int nzc20[7] = {20, 21, 25, 26, 27, 28, 29};
  
  double G[45];
  G[0] = 0.5*K[0]*K[1]*det;
  G[1] = 0.5*K[1]*K[3]*det;
  G[2] = 0.5*K[1]*K[6]*det;
  G[3] = det*(0.5*(K[0]*K[3] + K[2]*K[5]) + K[1]*K[4]);
  G[4] = 0.5*K[3]*K[4]*det;
  G[5] = det*(0.5*(K[3]*K[3] + K[5]*K[5]) + K[4]*K[4]);
  G[6] = 0.5*K[3]*K[7]*det;
  G[7] = det*(0.5*(K[3]*K[6] + K[5]*K[8]) + K[4]*K[7]);
  G[8] = det*(0.5*(K[0]*K[6] + K[2]*K[8]) + K[1]*K[7]);
  G[9] = 0.5*K[4]*K[6]*det;
  G[10] = 0.5*K[6]*K[7]*det;
  G[11] = det*(0.5*(K[6]*K[6] + K[8]*K[8]) + K[7]*K[7]);
  G[12] = det*(0.5*(K[1]*K[4] + K[2]*K[5]) + K[0]*K[3]);
  G[13] = 0.5*K[0]*K[4]*det;
  G[14] = det*(0.5*(K[0]*K[0] + K[2]*K[2]) + K[1]*K[1]);
  G[15] = 0.5*K[0]*K[7]*det;
  G[16] = det*(0.5*(K[4]*K[7] + K[5]*K[8]) + K[3]*K[6]);
  G[17] = det*(0.5*(K[1]*K[7] + K[2]*K[8]) + K[0]*K[6]);
  G[18] = det*(0.5*(K[7]*K[7] + K[8]*K[8]) + K[6]*K[6]);
  G[19] = det*(0.5*(K[4]*K[4] + K[5]*K[5]) + K[3]*K[3]);
  G[20] = 0.5*K[1]*K[5]*det;
  G[21] = 0.5*K[4]*K[5]*det;
  G[22] = 0.5*K[5]*K[7]*det;
  G[23] = 0.5*K[7]*K[8]*det;
  G[24] = 0.5*K[2]*K[7]*det;
  G[25] = 0.5*K[1]*K[8]*det;
  G[26] = 0.5*K[1]*K[2]*det;
  G[27] = 0.5*K[4]*K[8]*det;
  G[28] = 0.5*K[2]*K[4]*det;
  G[29] = det*(0.5*(K[1]*K[1] + K[2]*K[2]) + K[0]*K[0]);
  G[30] = 0.5*K[3]*K[5]*det;
  G[31] = 0.5*K[5]*K[6]*det;
  G[32] = 0.5*K[0]*K[5]*det;
  G[33] = 0.5*K[2]*K[3]*det;
  G[34] = 0.5*K[2]*K[6]*det;
  G[35] = 0.5*K[0]*K[2]*det;
  G[36] = 0.5*K[3]*K[8]*det;
  G[37] = 0.5*K[6]*K[8]*det;
  G[38] = 0.5*K[0]*K[8]*det;
  G[39] = det*(0.5*(K[3]*K[3] + K[4]*K[4]) + K[5]*K[5]);
  G[40] = det*(0.5*(K[3]*K[6] + K[4]*K[7]) + K[5]*K[8]);
  G[41] = det*(0.5*(K[0]*K[3] + K[1]*K[4]) + K[2]*K[5]);
  G[42] = det*(0.5*(K[0]*K[6] + K[1]*K[7]) + K[2]*K[8]);
  G[43] = det*(0.5*(K[0]*K[0] + K[1]*K[1]) + K[2]*K[2]);
  G[44] = det*(0.5*(K[6]*K[6] + K[7]*K[7]) + K[8]*K[8]);
  
  
  for (int ip  = 0; ip < 4; ip += 1)
  {
    double Z[45];
    Z[0] = G[0]*W4[ip];
    Z[1] = G[1]*W4[ip];
    Z[2] = G[2]*W4[ip];
    Z[3] = G[3]*W4[ip];
    Z[4] = G[4]*W4[ip];
    Z[5] = G[5]*W4[ip];
    Z[6] = G[6]*W4[ip];
    Z[7] = G[7]*W4[ip];
    Z[8] = G[8]*W4[ip];
    Z[9] = G[9]*W4[ip];
    Z[10] = G[10]*W4[ip];
    Z[11] = G[11]*W4[ip];
    Z[12] = G[12]*W4[ip];
    Z[13] = G[13]*W4[ip];
    Z[14] = G[14]*W4[ip];
    Z[15] = G[15]*W4[ip];
    Z[16] = G[16]*W4[ip];
    Z[17] = G[17]*W4[ip];
    Z[18] = G[18]*W4[ip];
    Z[19] = G[19]*W4[ip];
    Z[20] = G[20]*W4[ip];
    Z[21] = G[21]*W4[ip];
    Z[22] = G[22]*W4[ip];
    Z[23] = G[23]*W4[ip];
    Z[24] = G[24]*W4[ip];
    Z[25] = G[25]*W4[ip];
    Z[26] = G[26]*W4[ip];
    Z[27] = G[27]*W4[ip];
    Z[28] = G[28]*W4[ip];
    Z[29] = G[29]*W4[ip];
    Z[30] = G[30]*W4[ip];
    Z[31] = G[31]*W4[ip];
    Z[32] = G[32]*W4[ip];
    Z[33] = G[33]*W4[ip];
    Z[34] = G[34]*W4[ip];
    Z[35] = G[35]*W4[ip];
    Z[36] = G[36]*W4[ip];
    Z[37] = G[37]*W4[ip];
    Z[38] = G[38]*W4[ip];
    Z[39] = G[39]*W4[ip];
    Z[40] = G[40]*W4[ip];
    Z[41] = G[41]*W4[ip];
    Z[42] = G[42]*W4[ip];
    Z[43] = G[43]*W4[ip];
    Z[44] = G[44]*W4[ip];
    
    for (int j  = 0; j < 7; j += 1)
    {
      
      for (int k  = 0; k < 7; k += 1)
      {
        A[nzc14[j]][nzc14[k]] += (Z[11] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc14[j]][nzc15[k]] += (Z[7] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc14[j]][nzc16[k]] += (Z[8] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc14[j]][nzc18[k]] += (Z[23] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc14[j]][nzc19[k]] += (Z[27] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc14[j]][nzc20[k]] += (Z[25] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc14[j]][nzc21[k]] += (Z[10] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc14[j]][nzc22[k]] += (Z[9] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc14[j]][nzc23[k]] += (Z[2] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc15[j]][nzc14[k]] += (Z[7] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc15[j]][nzc15[k]] += (Z[5] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc15[j]][nzc16[k]] += (Z[3] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc15[j]][nzc18[k]] += (Z[22] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc15[j]][nzc19[k]] += (Z[21] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc15[j]][nzc20[k]] += (Z[20] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc15[j]][nzc21[k]] += (Z[6] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc15[j]][nzc22[k]] += (Z[4] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc15[j]][nzc23[k]] += (Z[1] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc16[j]][nzc14[k]] += (Z[8] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc16[j]][nzc15[k]] += (Z[3] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc16[j]][nzc16[k]] += (Z[14] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        A[nzc16[j]][nzc18[k]] += (Z[24] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc16[j]][nzc19[k]] += (Z[28] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc16[j]][nzc20[k]] += (Z[26] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        A[nzc16[j]][nzc21[k]] += (Z[15] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc16[j]][nzc22[k]] += (Z[13] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc16[j]][nzc23[k]] += (Z[0] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        A[nzc18[j]][nzc14[k]] += (Z[23] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc18[j]][nzc15[k]] += (Z[22] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc18[j]][nzc16[k]] += (Z[24] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc18[j]][nzc18[k]] += (Z[44] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc18[j]][nzc19[k]] += (Z[40] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc18[j]][nzc20[k]] += (Z[42] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc18[j]][nzc21[k]] += (Z[37] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc18[j]][nzc22[k]] += (Z[31] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc18[j]][nzc23[k]] += (Z[34] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc19[j]][nzc14[k]] += (Z[27] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc19[j]][nzc15[k]] += (Z[21] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc19[j]][nzc16[k]] += (Z[28] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc19[j]][nzc18[k]] += (Z[40] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc19[j]][nzc19[k]] += (Z[39] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc19[j]][nzc20[k]] += (Z[41] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc19[j]][nzc21[k]] += (Z[36] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc19[j]][nzc22[k]] += (Z[30] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc19[j]][nzc23[k]] += (Z[33] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc20[j]][nzc14[k]] += (Z[25] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc20[j]][nzc15[k]] += (Z[20] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc20[j]][nzc16[k]] += (Z[26] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        A[nzc20[j]][nzc18[k]] += (Z[42] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc20[j]][nzc19[k]] += (Z[41] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc20[j]][nzc20[k]] += (Z[43] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        A[nzc20[j]][nzc21[k]] += (Z[38] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc20[j]][nzc22[k]] += (Z[32] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc20[j]][nzc23[k]] += (Z[35] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        A[nzc21[j]][nzc14[k]] += (Z[10] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc21[j]][nzc15[k]] += (Z[6] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc21[j]][nzc16[k]] += (Z[15] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc21[j]][nzc18[k]] += (Z[37] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc21[j]][nzc19[k]] += (Z[36] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc21[j]][nzc20[k]] += (Z[38] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc21[j]][nzc21[k]] += (Z[18] * FE1_C1_D001[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc21[j]][nzc22[k]] += (Z[16] * FE1_C1_D010[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc21[j]][nzc23[k]] += (Z[17] * FE1_C1_D100[ip][k] * FE1_C1_D001[ip][j]);
        A[nzc22[j]][nzc14[k]] += (Z[9] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc22[j]][nzc15[k]] += (Z[4] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc22[j]][nzc16[k]] += (Z[13] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc22[j]][nzc18[k]] += (Z[31] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc22[j]][nzc19[k]] += (Z[30] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc22[j]][nzc20[k]] += (Z[32] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc22[j]][nzc21[k]] += (Z[16] * FE1_C1_D010[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc22[j]][nzc22[k]] += (Z[19] * FE1_C1_D010[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc22[j]][nzc23[k]] += (Z[12] * FE1_C1_D100[ip][k] * FE1_C1_D010[ip][j]);
        A[nzc23[j]][nzc14[k]] += (Z[2] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc23[j]][nzc15[k]] += (Z[1] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc23[j]][nzc16[k]] += (Z[0] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        A[nzc23[j]][nzc18[k]] += (Z[34] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc23[j]][nzc19[k]] += (Z[33] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc23[j]][nzc20[k]] += (Z[35] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        A[nzc23[j]][nzc21[k]] += (Z[17] * FE1_C1_D100[ip][j] * FE1_C1_D001[ip][k]);
        A[nzc23[j]][nzc22[k]] += (Z[12] * FE1_C1_D100[ip][j] * FE1_C1_D010[ip][k]);
        A[nzc23[j]][nzc23[k]] += (Z[29] * FE1_C1_D100[ip][k] * FE1_C1_D100[ip][j]);
        
      }
      
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
    double buffer_arg0_0[30][30]  = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec);
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
        
        