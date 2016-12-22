
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


void form_cell_integral_0_otherwise(double A[4][4], double *coordinate_dofs)
{
// Section for quadrature weights and points
static const double weights5[5] = { -0.133333333333, 0.075, 0.075, 0.075, 0.075 };
static const double points5[15] = { 0.25, 0.25, 0.25, 0.5, 0.166666666667, 0.166666666667, 0.166666666667, 0.5, 0.166666666667, 0.166666666667, 0.166666666667, 0.5, 0.166666666667, 0.166666666667, 0.166666666667 };
// Section for precomputed element basis function values
// Table dimensions: num_entities, num_points, num_dofs
// Definitions of 4 tables for 5 quadrature points
static const double FE0_C0_D001_Q5[1][5][4] =
    { { { -1.0, -5.55111512313e-17, 0.0, 1.0 },
        { -1.0, -1.11022302463e-16, 0.0, 1.0 },
        { -1.0, -1.11022302463e-16, 0.0, 1.0 },
        { -1.0, -1.66533453694e-16, -5.55111512313e-17, 1.0 },
        { -1.0, -1.11022302463e-16, 0.0, 1.0 } } };
static const double FE0_C0_D010_Q5[1][5][3] =
    { { { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -1.66533453694e-16, 1.0 },
        { -1.0, -1.66533453694e-16, 1.0 },
        { -1.0, -1.66533453694e-16, 1.0 },
        { -1.0, -1.66533453694e-16, 1.0 } } };
static const double FE0_C0_D100_Q5[1][5][2] =
    { { { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 } } };
static const double FE1_C0_Q5[1][5][4] =
    { { { 0.25, 0.25, 0.25, 0.25 },
        { 0.166666666667, 0.5, 0.166666666667, 0.166666666667 },
        { 0.166666666667, 0.166666666667, 0.5, 0.166666666667 },
        { 0.166666666667, 0.166666666667, 0.166666666667, 0.5 },
        { 0.5, 0.166666666667, 0.166666666667, 0.166666666667 } } };
// Section for piecewise constant computations
const double J_c4 = coordinate_dofs[1] * FE0_C0_D010_Q5[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D010_Q5[0][0][5 - 4] + coordinate_dofs[7] * FE0_C0_D010_Q5[0][0][6 - 4];
const double J_c8 = coordinate_dofs[2] * FE0_C0_D001_Q5[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D001_Q5[0][0][9 - 8] + coordinate_dofs[8] * FE0_C0_D001_Q5[0][0][10 - 8] + coordinate_dofs[11] * FE0_C0_D001_Q5[0][0][11 - 8];
const double J_c5 = coordinate_dofs[1] * FE0_C0_D001_Q5[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D001_Q5[0][0][5 - 4] + coordinate_dofs[7] * FE0_C0_D001_Q5[0][0][6 - 4] + coordinate_dofs[10] * FE0_C0_D001_Q5[0][0][7 - 4];
const double J_c7 = coordinate_dofs[2] * FE0_C0_D010_Q5[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D010_Q5[0][0][9 - 8] + coordinate_dofs[8] * FE0_C0_D010_Q5[0][0][10 - 8];
const double J_c0 = coordinate_dofs[0] * FE0_C0_D100_Q5[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D100_Q5[0][0][1 - 0];
const double J_c1 = coordinate_dofs[0] * FE0_C0_D010_Q5[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D010_Q5[0][0][1 - 0] + coordinate_dofs[6] * FE0_C0_D010_Q5[0][0][2 - 0];
const double J_c6 = coordinate_dofs[2] * FE0_C0_D100_Q5[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D100_Q5[0][0][9 - 8];
const double J_c3 = coordinate_dofs[1] * FE0_C0_D100_Q5[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D100_Q5[0][0][5 - 4];
const double J_c2 = coordinate_dofs[0] * FE0_C0_D001_Q5[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D001_Q5[0][0][1 - 0] + coordinate_dofs[6] * FE0_C0_D001_Q5[0][0][2 - 0] + coordinate_dofs[9] * FE0_C0_D001_Q5[0][0][3 - 0];
double sp5[84];
sp5[0] = J_c4 * J_c8;
sp5[1] = J_c5 * J_c7;
sp5[2] = -1 * sp5[1];
sp5[3] = sp5[0] + sp5[2];
sp5[4] = J_c0 * sp5[3];
sp5[5] = J_c5 * J_c6;
sp5[6] = J_c3 * J_c8;
sp5[7] = -1 * sp5[6];
sp5[8] = sp5[5] + sp5[7];
sp5[9] = J_c1 * sp5[8];
sp5[10] = sp5[4] + sp5[9];
sp5[11] = J_c3 * J_c7;
sp5[12] = J_c4 * J_c6;
sp5[13] = -1 * sp5[12];
sp5[14] = sp5[11] + sp5[13];
sp5[15] = J_c2 * sp5[14];
sp5[16] = sp5[10] + sp5[15];
sp5[17] = sp5[3] / sp5[16];
sp5[18] = -1 * J_c8;
sp5[19] = J_c3 * sp5[18];
sp5[20] = sp5[5] + sp5[19];
sp5[21] = sp5[20] / sp5[16];
sp5[22] = sp5[14] / sp5[16];
sp5[23] = sp5[22] * sp5[21];
sp5[24] = sp5[22] * sp5[22];
sp5[25] = sp5[22] * sp5[17];
sp5[26] = sp5[17] * sp5[21];
sp5[27] = sp5[17] * sp5[17];
sp5[28] = sp5[21] * sp5[21];
sp5[29] = J_c2 * J_c7;
sp5[30] = -1 * J_c1;
sp5[31] = J_c8 * sp5[30];
sp5[32] = sp5[29] + sp5[31];
sp5[33] = sp5[32] / sp5[16];
sp5[34] = J_c0 * J_c8;
sp5[35] = -1 * J_c2;
sp5[36] = J_c6 * sp5[35];
sp5[37] = sp5[34] + sp5[36];
sp5[38] = sp5[37] / sp5[16];
sp5[39] = J_c1 * J_c6;
sp5[40] = J_c0 * J_c7;
sp5[41] = -1 * sp5[40];
sp5[42] = sp5[39] + sp5[41];
sp5[43] = sp5[42] / sp5[16];
sp5[44] = sp5[38] * sp5[43];
sp5[45] = sp5[43] * sp5[43];
sp5[46] = sp5[33] * sp5[43];
sp5[47] = sp5[33] * sp5[38];
sp5[48] = sp5[33] * sp5[33];
sp5[49] = sp5[38] * sp5[38];
sp5[50] = sp5[48] + sp5[27];
sp5[51] = sp5[47] + sp5[26];
sp5[52] = sp5[46] + sp5[25];
sp5[53] = sp5[49] + sp5[28];
sp5[54] = sp5[44] + sp5[23];
sp5[55] = sp5[24] + sp5[45];
sp5[56] = J_c1 * J_c5;
sp5[57] = J_c2 * J_c4;
sp5[58] = -1 * sp5[57];
sp5[59] = sp5[56] + sp5[58];
sp5[60] = sp5[59] / sp5[16];
sp5[61] = J_c2 * J_c3;
sp5[62] = J_c0 * J_c5;
sp5[63] = -1 * sp5[62];
sp5[64] = sp5[61] + sp5[63];
sp5[65] = sp5[64] / sp5[16];
sp5[66] = J_c0 * J_c4;
sp5[67] = J_c1 * J_c3;
sp5[68] = -1 * sp5[67];
sp5[69] = sp5[66] + sp5[68];
sp5[70] = sp5[69] / sp5[16];
sp5[71] = sp5[70] * sp5[65];
sp5[72] = sp5[70] * sp5[70];
sp5[73] = sp5[70] * sp5[60];
sp5[74] = sp5[60] * sp5[65];
sp5[75] = sp5[60] * sp5[60];
sp5[76] = sp5[65] * sp5[65];
sp5[77] = sp5[50] + sp5[75];
sp5[78] = sp5[51] + sp5[74];
sp5[79] = sp5[52] + sp5[73];
sp5[80] = sp5[53] + sp5[76];
sp5[81] = sp5[54] + sp5[71];
sp5[82] = sp5[55] + sp5[72];
sp5[83] = fabs(sp5[16]);
for (int iq = 0; iq < 5; ++iq)
{
    // Quadrature loop body setup (num_points=5)
    // Section for geometrically varying computations
    double sv5[7];
    sv5[0] = weights5[iq] * sp5[83];
    sv5[1] = sp5[79] * sv5[0];
    sv5[2] = sp5[81] * sv5[0];
    sv5[3] = sp5[78] * sv5[0];
    sv5[4] = sp5[80] * sv5[0];
    sv5[5] = sp5[77] * sv5[0];
    sv5[6] = sp5[82] * sv5[0];
    for (int ia0 = 0; ia0 < 2; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[5] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[3] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[1] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
    }
    for (int ia0 = 0; ia0 < 3; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[3] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[4] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[2] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
    }
    for (int ia0 = 0; ia0 < 4; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[1] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[2] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[0] * FE1_C0_Q5[0][iq][ia0 - 0] * FE1_C0_Q5[0][iq][ia1 - 0];
            A[ia0][ia1] += sv5[6] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
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
    double buffer_arg0_0[4][4] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec);
    MatSetValuesLocal(arg0_0_0, 4, arg0_0_map0_0 + i * 4,
                                             4, arg0_0_map1_0 + i * 4,
                                             (const PetscScalar *)buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        