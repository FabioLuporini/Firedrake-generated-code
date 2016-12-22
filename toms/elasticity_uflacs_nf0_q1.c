
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


void form_cell_integral_0_otherwise(double A[12][12], double *coordinate_dofs)
{
// Section for quadrature weights and points
static const double weights1[1] = { 0.166666666667 };
static const double points1[3] = { 0.25, 0.25, 0.25 };
// Section for precomputed element basis function values
// Table dimensions: num_entities, num_points, num_dofs
// Definitions of 3 tables for 1 quadrature points
static const double FE0_C0_D001_Q1[1][1][4] = { { { -1.0, -5.55111512313e-17, 0.0, 1.0 } } };
static const double FE0_C0_D010_Q1[1][1][3] = { { { -1.0, -1.11022302463e-16, 1.0 } } };
static const double FE0_C0_D100_Q1[1][1][2] = { { { -1.0, 1.0 } } };
// Section for piecewise constant computations
const double J_c4 = coordinate_dofs[1] * FE0_C0_D010_Q1[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D010_Q1[0][0][5 - 4] + coordinate_dofs[7] * FE0_C0_D010_Q1[0][0][6 - 4];
const double J_c8 = coordinate_dofs[2] * FE0_C0_D001_Q1[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D001_Q1[0][0][9 - 8] + coordinate_dofs[8] * FE0_C0_D001_Q1[0][0][10 - 8] + coordinate_dofs[11] * FE0_C0_D001_Q1[0][0][11 - 8];
const double J_c5 = coordinate_dofs[1] * FE0_C0_D001_Q1[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D001_Q1[0][0][5 - 4] + coordinate_dofs[7] * FE0_C0_D001_Q1[0][0][6 - 4] + coordinate_dofs[10] * FE0_C0_D001_Q1[0][0][7 - 4];
const double J_c7 = coordinate_dofs[2] * FE0_C0_D010_Q1[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D010_Q1[0][0][9 - 8] + coordinate_dofs[8] * FE0_C0_D010_Q1[0][0][10 - 8];
const double J_c0 = coordinate_dofs[0] * FE0_C0_D100_Q1[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D100_Q1[0][0][1 - 0];
const double J_c1 = coordinate_dofs[0] * FE0_C0_D010_Q1[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D010_Q1[0][0][1 - 0] + coordinate_dofs[6] * FE0_C0_D010_Q1[0][0][2 - 0];
const double J_c6 = coordinate_dofs[2] * FE0_C0_D100_Q1[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D100_Q1[0][0][9 - 8];
const double J_c3 = coordinate_dofs[1] * FE0_C0_D100_Q1[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D100_Q1[0][0][5 - 4];
const double J_c2 = coordinate_dofs[0] * FE0_C0_D001_Q1[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D001_Q1[0][0][1 - 0] + coordinate_dofs[6] * FE0_C0_D001_Q1[0][0][2 - 0] + coordinate_dofs[9] * FE0_C0_D001_Q1[0][0][3 - 0];
double sp1[264];
sp1[0] = J_c4 * J_c8;
sp1[1] = J_c5 * J_c7;
sp1[2] = -1 * sp1[1];
sp1[3] = sp1[0] + sp1[2];
sp1[4] = J_c0 * sp1[3];
sp1[5] = J_c5 * J_c6;
sp1[6] = J_c3 * J_c8;
sp1[7] = -1 * sp1[6];
sp1[8] = sp1[5] + sp1[7];
sp1[9] = J_c1 * sp1[8];
sp1[10] = sp1[4] + sp1[9];
sp1[11] = J_c3 * J_c7;
sp1[12] = J_c4 * J_c6;
sp1[13] = -1 * sp1[12];
sp1[14] = sp1[11] + sp1[13];
sp1[15] = J_c2 * sp1[14];
sp1[16] = sp1[10] + sp1[15];
sp1[17] = sp1[3] / sp1[16];
sp1[18] = -1 * J_c8;
sp1[19] = J_c3 * sp1[18];
sp1[20] = sp1[5] + sp1[19];
sp1[21] = sp1[20] / sp1[16];
sp1[22] = sp1[14] / sp1[16];
sp1[23] = sp1[17] + sp1[17];
sp1[24] = sp1[21] + sp1[21];
sp1[25] = sp1[22] + sp1[22];
sp1[26] = sp1[24] * sp1[24];
sp1[27] = sp1[25] * sp1[24];
sp1[28] = sp1[23] * sp1[24];
sp1[29] = sp1[25] * sp1[25];
sp1[30] = sp1[25] * sp1[23];
sp1[31] = sp1[23] * sp1[23];
sp1[32] = J_c2 * J_c7;
sp1[33] = -1 * J_c1;
sp1[34] = J_c8 * sp1[33];
sp1[35] = sp1[32] + sp1[34];
sp1[36] = sp1[35] / sp1[16];
sp1[37] = J_c0 * J_c8;
sp1[38] = -1 * J_c2;
sp1[39] = J_c6 * sp1[38];
sp1[40] = sp1[37] + sp1[39];
sp1[41] = sp1[40] / sp1[16];
sp1[42] = J_c1 * J_c6;
sp1[43] = J_c0 * J_c7;
sp1[44] = -1 * sp1[43];
sp1[45] = sp1[42] + sp1[44];
sp1[46] = sp1[45] / sp1[16];
sp1[47] = sp1[22] * sp1[22];
sp1[48] = sp1[22] * sp1[46];
sp1[49] = sp1[36] * sp1[22];
sp1[50] = sp1[22] * sp1[17];
sp1[51] = sp1[41] * sp1[22];
sp1[52] = sp1[22] * sp1[21];
sp1[53] = sp1[46] * sp1[17];
sp1[54] = sp1[36] * sp1[17];
sp1[55] = sp1[17] * sp1[17];
sp1[56] = sp1[41] * sp1[17];
sp1[57] = sp1[17] * sp1[21];
sp1[58] = sp1[36] * sp1[46];
sp1[59] = sp1[36] * sp1[36];
sp1[60] = sp1[36] * sp1[41];
sp1[61] = sp1[36] * sp1[21];
sp1[62] = sp1[46] * sp1[21];
sp1[63] = sp1[41] * sp1[21];
sp1[64] = sp1[21] * sp1[21];
sp1[65] = sp1[41] * sp1[46];
sp1[66] = sp1[41] * sp1[41];
sp1[67] = sp1[46] * sp1[46];
sp1[68] = sp1[31] + sp1[59];
sp1[69] = sp1[28] + sp1[60];
sp1[70] = sp1[30] + sp1[58];
sp1[71] = sp1[26] + sp1[66];
sp1[72] = sp1[27] + sp1[65];
sp1[73] = sp1[29] + sp1[67];
sp1[74] = J_c1 * J_c5;
sp1[75] = J_c2 * J_c4;
sp1[76] = -1 * sp1[75];
sp1[77] = sp1[74] + sp1[76];
sp1[78] = sp1[77] / sp1[16];
sp1[79] = J_c2 * J_c3;
sp1[80] = J_c0 * J_c5;
sp1[81] = -1 * sp1[80];
sp1[82] = sp1[79] + sp1[81];
sp1[83] = sp1[82] / sp1[16];
sp1[84] = J_c0 * J_c4;
sp1[85] = J_c1 * J_c3;
sp1[86] = -1 * sp1[85];
sp1[87] = sp1[84] + sp1[86];
sp1[88] = sp1[87] / sp1[16];
sp1[89] = sp1[88] * sp1[22];
sp1[90] = sp1[78] * sp1[22];
sp1[91] = sp1[83] * sp1[22];
sp1[92] = sp1[88] * sp1[17];
sp1[93] = sp1[78] * sp1[17];
sp1[94] = sp1[83] * sp1[17];
sp1[95] = sp1[88] * sp1[78];
sp1[96] = sp1[78] * sp1[78];
sp1[97] = sp1[78] * sp1[21];
sp1[98] = sp1[78] * sp1[83];
sp1[99] = sp1[88] * sp1[83];
sp1[100] = sp1[83] * sp1[21];
sp1[101] = sp1[83] * sp1[83];
sp1[102] = sp1[88] * sp1[88];
sp1[103] = sp1[88] * sp1[21];
sp1[104] = sp1[68] + sp1[96];
sp1[105] = sp1[69] + sp1[98];
sp1[106] = sp1[70] + sp1[95];
sp1[107] = sp1[71] + sp1[101];
sp1[108] = sp1[72] + sp1[99];
sp1[109] = sp1[73] + sp1[102];
sp1[110] = sp1[36] + sp1[36];
sp1[111] = sp1[41] + sp1[41];
sp1[112] = sp1[46] + sp1[46];
sp1[113] = sp1[111] * sp1[112];
sp1[114] = sp1[110] * sp1[111];
sp1[115] = sp1[111] * sp1[111];
sp1[116] = sp1[112] * sp1[112];
sp1[117] = sp1[110] * sp1[112];
sp1[118] = sp1[110] * sp1[110];
sp1[119] = sp1[118] + sp1[55];
sp1[120] = sp1[114] + sp1[57];
sp1[121] = sp1[117] + sp1[50];
sp1[122] = sp1[115] + sp1[64];
sp1[123] = sp1[113] + sp1[52];
sp1[124] = sp1[116] + sp1[47];
sp1[125] = sp1[88] * sp1[46];
sp1[126] = sp1[36] * sp1[88];
sp1[127] = sp1[41] * sp1[88];
sp1[128] = sp1[78] * sp1[46];
sp1[129] = sp1[36] * sp1[78];
sp1[130] = sp1[41] * sp1[78];
sp1[131] = sp1[83] * sp1[46];
sp1[132] = sp1[36] * sp1[83];
sp1[133] = sp1[41] * sp1[83];
sp1[134] = sp1[119] + sp1[96];
sp1[135] = sp1[120] + sp1[98];
sp1[136] = sp1[121] + sp1[95];
sp1[137] = sp1[122] + sp1[101];
sp1[138] = sp1[123] + sp1[99];
sp1[139] = sp1[124] + sp1[102];
sp1[140] = sp1[104] + sp1[59];
sp1[141] = sp1[105] + sp1[60];
sp1[142] = sp1[106] + sp1[58];
sp1[143] = sp1[54] + sp1[54];
sp1[144] = sp1[61] + sp1[61];
sp1[145] = sp1[49] + sp1[49];
sp1[146] = sp1[107] + sp1[66];
sp1[147] = sp1[108] + sp1[65];
sp1[148] = sp1[56] + sp1[56];
sp1[149] = sp1[63] + sp1[63];
sp1[150] = sp1[51] + sp1[51];
sp1[151] = sp1[109] + sp1[67];
sp1[152] = sp1[53] + sp1[53];
sp1[153] = sp1[62] + sp1[62];
sp1[154] = sp1[48] + sp1[48];
sp1[155] = sp1[134] + sp1[55];
sp1[156] = sp1[135] + sp1[57];
sp1[157] = sp1[136] + sp1[50];
sp1[158] = sp1[137] + sp1[64];
sp1[159] = sp1[138] + sp1[52];
sp1[160] = sp1[139] + sp1[47];
sp1[161] = sp1[59] + sp1[55];
sp1[162] = sp1[60] + sp1[57];
sp1[163] = sp1[58] + sp1[50];
sp1[164] = sp1[66] + sp1[64];
sp1[165] = sp1[65] + sp1[52];
sp1[166] = sp1[47] + sp1[67];
sp1[167] = sp1[78] + sp1[78];
sp1[168] = sp1[83] + sp1[83];
sp1[169] = sp1[88] + sp1[88];
sp1[170] = sp1[168] * sp1[168];
sp1[171] = sp1[169] * sp1[168];
sp1[172] = sp1[167] * sp1[168];
sp1[173] = sp1[169] * sp1[169];
sp1[174] = sp1[169] * sp1[167];
sp1[175] = sp1[167] * sp1[167];
sp1[176] = sp1[161] + sp1[175];
sp1[177] = sp1[162] + sp1[172];
sp1[178] = sp1[163] + sp1[174];
sp1[179] = sp1[164] + sp1[170];
sp1[180] = sp1[165] + sp1[171];
sp1[181] = sp1[166] + sp1[173];
sp1[182] = sp1[140] + sp1[96];
sp1[183] = sp1[141] + sp1[98];
sp1[184] = sp1[142] + sp1[95];
sp1[185] = sp1[93] + sp1[93];
sp1[186] = sp1[97] + sp1[97];
sp1[187] = sp1[90] + sp1[90];
sp1[188] = sp1[146] + sp1[101];
sp1[189] = sp1[147] + sp1[99];
sp1[190] = sp1[94] + sp1[94];
sp1[191] = sp1[100] + sp1[100];
sp1[192] = sp1[91] + sp1[91];
sp1[193] = sp1[151] + sp1[102];
sp1[194] = sp1[92] + sp1[92];
sp1[195] = sp1[103] + sp1[103];
sp1[196] = sp1[89] + sp1[89];
sp1[197] = sp1[155] + sp1[96];
sp1[198] = sp1[156] + sp1[98];
sp1[199] = sp1[157] + sp1[95];
sp1[200] = sp1[129] + sp1[129];
sp1[201] = sp1[130] + sp1[130];
sp1[202] = sp1[128] + sp1[128];
sp1[203] = sp1[158] + sp1[101];
sp1[204] = sp1[159] + sp1[99];
sp1[205] = sp1[132] + sp1[132];
sp1[206] = sp1[133] + sp1[133];
sp1[207] = sp1[131] + sp1[131];
sp1[208] = sp1[160] + sp1[102];
sp1[209] = sp1[126] + sp1[126];
sp1[210] = sp1[127] + sp1[127];
sp1[211] = sp1[125] + sp1[125];
sp1[212] = sp1[176] + sp1[161];
sp1[213] = sp1[177] + sp1[162];
sp1[214] = sp1[178] + sp1[163];
sp1[215] = sp1[179] + sp1[164];
sp1[216] = sp1[180] + sp1[165];
sp1[217] = sp1[181] + sp1[166];
sp1[218] = 0.25 * sp1[216];
sp1[219] = 0.25 * sp1[211];
sp1[220] = 0.25 * sp1[191];
sp1[221] = 0.25 * sp1[199];
sp1[222] = 0.25 * sp1[204];
sp1[223] = 0.25 * sp1[145];
sp1[224] = 0.25 * sp1[184];
sp1[225] = 0.25 * sp1[150];
sp1[226] = 0.25 * sp1[198];
sp1[227] = 0.25 * sp1[153];
sp1[228] = 0.25 * sp1[207];
sp1[229] = 0.25 * sp1[144];
sp1[230] = 0.25 * sp1[203];
sp1[231] = 0.25 * sp1[186];
sp1[232] = 0.25 * sp1[154];
sp1[233] = 0.25 * sp1[195];
sp1[234] = 0.25 * sp1[206];
sp1[235] = 0.25 * sp1[148];
sp1[236] = 0.25 * sp1[205];
sp1[237] = 0.25 * sp1[183];
sp1[238] = 0.25 * sp1[196];
sp1[239] = 0.25 * sp1[152];
sp1[240] = 0.25 * sp1[209];
sp1[241] = 0.25 * sp1[214];
sp1[242] = 0.25 * sp1[192];
sp1[243] = 0.25 * sp1[200];
sp1[244] = 0.25 * sp1[143];
sp1[245] = 0.25 * sp1[194];
sp1[246] = 0.25 * sp1[201];
sp1[247] = 0.25 * sp1[189];
sp1[248] = 0.25 * sp1[149];
sp1[249] = 0.25 * sp1[210];
sp1[250] = 0.25 * sp1[197];
sp1[251] = 0.25 * sp1[193];
sp1[252] = 0.25 * sp1[190];
sp1[253] = 0.25 * sp1[202];
sp1[254] = 0.25 * sp1[213];
sp1[255] = 0.25 * sp1[185];
sp1[256] = 0.25 * sp1[187];
sp1[257] = 0.25 * sp1[217];
sp1[258] = 0.25 * sp1[215];
sp1[259] = 0.25 * sp1[208];
sp1[260] = 0.25 * sp1[212];
sp1[261] = 0.25 * sp1[188];
sp1[262] = 0.25 * sp1[182];
sp1[263] = fabs(sp1[16]);
// Only 1 quadrature point, no loop
const int iq = 0;
{
    // Quadrature loop body setup (num_points=1)
    // Section for geometrically varying computations
    double sv1[46];
    sv1[0] = weights1[iq] * sp1[263];
    sv1[1] = sp1[218] * sv1[0];
    sv1[2] = sp1[220] * sv1[0];
    sv1[3] = sp1[257] * sv1[0];
    sv1[4] = sp1[222] * sv1[0];
    sv1[5] = sp1[223] * sv1[0];
    sv1[6] = sp1[224] * sv1[0];
    sv1[7] = sp1[246] * sv1[0];
    sv1[8] = sp1[225] * sv1[0];
    sv1[9] = sp1[226] * sv1[0];
    sv1[10] = sp1[227] * sv1[0];
    sv1[11] = sp1[228] * sv1[0];
    sv1[12] = sp1[229] * sv1[0];
    sv1[13] = sp1[231] * sv1[0];
    sv1[14] = sp1[232] * sv1[0];
    sv1[15] = sp1[219] * sv1[0];
    sv1[16] = sp1[233] * sv1[0];
    sv1[17] = sp1[234] * sv1[0];
    sv1[18] = sp1[235] * sv1[0];
    sv1[19] = sp1[236] * sv1[0];
    sv1[20] = sp1[237] * sv1[0];
    sv1[21] = sp1[238] * sv1[0];
    sv1[22] = sp1[239] * sv1[0];
    sv1[23] = sp1[240] * sv1[0];
    sv1[24] = sp1[221] * sv1[0];
    sv1[25] = sp1[243] * sv1[0];
    sv1[26] = sp1[242] * sv1[0];
    sv1[27] = sp1[230] * sv1[0];
    sv1[28] = sp1[245] * sv1[0];
    sv1[29] = sp1[256] * sv1[0];
    sv1[30] = sp1[247] * sv1[0];
    sv1[31] = sp1[241] * sv1[0];
    sv1[32] = sp1[249] * sv1[0];
    sv1[33] = sp1[250] * sv1[0];
    sv1[34] = sp1[251] * sv1[0];
    sv1[35] = sp1[252] * sv1[0];
    sv1[36] = sp1[253] * sv1[0];
    sv1[37] = sp1[255] * sv1[0];
    sv1[38] = sp1[244] * sv1[0];
    sv1[39] = sp1[248] * sv1[0];
    sv1[40] = sp1[258] * sv1[0];
    sv1[41] = sp1[259] * sv1[0];
    sv1[42] = sp1[260] * sv1[0];
    sv1[43] = sp1[254] * sv1[0];
    sv1[44] = sp1[261] * sv1[0];
    sv1[45] = sp1[262] * sv1[0];
    for (int ia0 = 0; ia0 < 2; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[45] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[20] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[6] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[38] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[12] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[5] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[37] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[13] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[29] * FE0_C0_D100_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 0; ia0 < 3; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[20] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[44] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[30] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[18] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[39] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[8] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[35] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[2] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[26] * FE0_C0_D010_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 0; ia0 < 4; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[6] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[30] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[34] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[22] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[10] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[14] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[28] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[16] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[21] * FE0_C0_D001_Q1[0][iq][ia0 - 0] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 6; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[38] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[18] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[22] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[33] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[9] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[24] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[25] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[7] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[36] * FE0_C0_D100_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 7; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[12] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[39] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[10] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[9] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[27] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[4] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[19] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[17] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[11] * FE0_C0_D010_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 8; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[5] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[8] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[14] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[24] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[4] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[41] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[23] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[32] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[15] * FE0_C0_D001_Q1[0][iq][ia0 - 4] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 10; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[37] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[35] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[28] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[25] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[19] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[23] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[42] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[43] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[31] * FE0_C0_D100_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 11; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[13] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[2] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[16] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[7] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[17] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[32] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[43] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[40] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[1] * FE0_C0_D010_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 12; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv1[29] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv1[26] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv1[21] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv1[36] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv1[11] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv1[15] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv1[31] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D100_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv1[1] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D010_Q1[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv1[3] * FE0_C0_D001_Q1[0][iq][ia0 - 8] * FE0_C0_D001_Q1[0][iq][ia1 - 8];
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
    double buffer_arg0_0[12][12] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec);
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
        
        