
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


void form_cell_integral_0_otherwise(double A[12][12], double *coordinate_dofs, double **w0)
{
// Section for quadrature weights and points
static const double weights4[4] = { 0.0416666666667, 0.0416666666667, 0.0416666666667, 0.0416666666667 };
static const double points4[12] = { 0.585410196625, 0.138196601125, 0.138196601125, 0.138196601125, 0.585410196625, 0.138196601125, 0.138196601125, 0.138196601125, 0.585410196625, 0.138196601125, 0.138196601125, 0.138196601125 };
// Section for precomputed element basis function values
// Table dimensions: num_entities, num_points, num_dofs
// Definitions of 4 tables for 4 quadrature points
static const double FE0_C0_D001_Q4[1][4][4] =
    { { { -1.0, -5.55111512313e-17, 0.0, 1.0 },
        { -1.0, -5.55111512313e-17, 0.0, 1.0 },
        { -1.0, -3.33066907388e-16, -1.11022302463e-16, 1.0 },
        { -1.0, -5.55111512313e-17, 0.0, 1.0 } } };
static const double FE0_C0_D010_Q4[1][4][3] =
    { { { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -2.77555756156e-16, 1.0 },
        { -1.0, -1.11022302463e-16, 1.0 } } };
static const double FE0_C0_D100_Q4[1][4][2] =
    { { { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 } } };
static const double FE1_C0_Q4[1][4][4] =
    { { { 0.138196601125, 0.585410196625, 0.138196601125, 0.138196601125 },
        { 0.138196601125, 0.138196601125, 0.585410196625, 0.138196601125 },
        { 0.138196601125, 0.138196601125, 0.138196601125, 0.585410196625 },
        { 0.585410196625, 0.138196601125, 0.138196601125, 0.138196601125 } } };
// Section for piecewise constant computations
const double J_c0 = coordinate_dofs[0] * FE0_C0_D100_Q4[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D100_Q4[0][0][1 - 0];
const double J_c4 = coordinate_dofs[1] * FE0_C0_D010_Q4[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D010_Q4[0][0][5 - 4] + coordinate_dofs[7] * FE0_C0_D010_Q4[0][0][6 - 4];
const double J_c8 = coordinate_dofs[2] * FE0_C0_D001_Q4[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D001_Q4[0][0][9 - 8] + coordinate_dofs[8] * FE0_C0_D001_Q4[0][0][10 - 8] + coordinate_dofs[11] * FE0_C0_D001_Q4[0][0][11 - 8];
const double J_c5 = coordinate_dofs[1] * FE0_C0_D001_Q4[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D001_Q4[0][0][5 - 4] + coordinate_dofs[7] * FE0_C0_D001_Q4[0][0][6 - 4] + coordinate_dofs[10] * FE0_C0_D001_Q4[0][0][7 - 4];
const double J_c7 = coordinate_dofs[2] * FE0_C0_D010_Q4[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D010_Q4[0][0][9 - 8] + coordinate_dofs[8] * FE0_C0_D010_Q4[0][0][10 - 8];
const double J_c1 = coordinate_dofs[0] * FE0_C0_D010_Q4[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D010_Q4[0][0][1 - 0] + coordinate_dofs[6] * FE0_C0_D010_Q4[0][0][2 - 0];
const double J_c6 = coordinate_dofs[2] * FE0_C0_D100_Q4[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D100_Q4[0][0][9 - 8];
const double J_c3 = coordinate_dofs[1] * FE0_C0_D100_Q4[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D100_Q4[0][0][5 - 4];
const double J_c2 = coordinate_dofs[0] * FE0_C0_D001_Q4[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D001_Q4[0][0][1 - 0] + coordinate_dofs[6] * FE0_C0_D001_Q4[0][0][2 - 0] + coordinate_dofs[9] * FE0_C0_D001_Q4[0][0][3 - 0];
double sp4[264];
sp4[0] = J_c4 * J_c8;
sp4[1] = J_c5 * J_c7;
sp4[2] = -1 * sp4[1];
sp4[3] = sp4[0] + sp4[2];
sp4[4] = J_c0 * sp4[3];
sp4[5] = J_c5 * J_c6;
sp4[6] = J_c3 * J_c8;
sp4[7] = -1 * sp4[6];
sp4[8] = sp4[5] + sp4[7];
sp4[9] = J_c1 * sp4[8];
sp4[10] = sp4[4] + sp4[9];
sp4[11] = J_c3 * J_c7;
sp4[12] = J_c4 * J_c6;
sp4[13] = -1 * sp4[12];
sp4[14] = sp4[11] + sp4[13];
sp4[15] = J_c2 * sp4[14];
sp4[16] = sp4[10] + sp4[15];
sp4[17] = fabs(sp4[16]);
sp4[18] = sp4[3] / sp4[16];
sp4[19] = -1 * J_c8;
sp4[20] = J_c3 * sp4[19];
sp4[21] = sp4[5] + sp4[20];
sp4[22] = sp4[21] / sp4[16];
sp4[23] = sp4[14] / sp4[16];
sp4[24] = sp4[18] + sp4[18];
sp4[25] = sp4[22] + sp4[22];
sp4[26] = sp4[23] + sp4[23];
sp4[27] = sp4[24] * sp4[25];
sp4[28] = sp4[26] * sp4[24];
sp4[29] = sp4[24] * sp4[24];
sp4[30] = sp4[25] * sp4[25];
sp4[31] = sp4[26] * sp4[25];
sp4[32] = sp4[26] * sp4[26];
sp4[33] = J_c2 * J_c7;
sp4[34] = -1 * J_c1;
sp4[35] = J_c8 * sp4[34];
sp4[36] = sp4[33] + sp4[35];
sp4[37] = sp4[36] / sp4[16];
sp4[38] = J_c0 * J_c8;
sp4[39] = -1 * J_c2;
sp4[40] = J_c6 * sp4[39];
sp4[41] = sp4[38] + sp4[40];
sp4[42] = sp4[41] / sp4[16];
sp4[43] = J_c1 * J_c6;
sp4[44] = J_c0 * J_c7;
sp4[45] = -1 * sp4[44];
sp4[46] = sp4[43] + sp4[45];
sp4[47] = sp4[46] / sp4[16];
sp4[48] = sp4[23] * sp4[22];
sp4[49] = sp4[18] * sp4[22];
sp4[50] = sp4[22] * sp4[22];
sp4[51] = sp4[42] * sp4[22];
sp4[52] = sp4[47] * sp4[22];
sp4[53] = sp4[37] * sp4[22];
sp4[54] = sp4[42] * sp4[23];
sp4[55] = sp4[42] * sp4[18];
sp4[56] = sp4[42] * sp4[42];
sp4[57] = sp4[42] * sp4[47];
sp4[58] = sp4[37] * sp4[42];
sp4[59] = sp4[23] * sp4[23];
sp4[60] = sp4[23] * sp4[18];
sp4[61] = sp4[23] * sp4[47];
sp4[62] = sp4[37] * sp4[23];
sp4[63] = sp4[47] * sp4[18];
sp4[64] = sp4[47] * sp4[47];
sp4[65] = sp4[37] * sp4[47];
sp4[66] = sp4[37] * sp4[18];
sp4[67] = sp4[37] * sp4[37];
sp4[68] = sp4[18] * sp4[18];
sp4[69] = sp4[29] + sp4[67];
sp4[70] = sp4[27] + sp4[58];
sp4[71] = sp4[28] + sp4[65];
sp4[72] = sp4[30] + sp4[56];
sp4[73] = sp4[31] + sp4[57];
sp4[74] = sp4[32] + sp4[64];
sp4[75] = J_c1 * J_c5;
sp4[76] = J_c2 * J_c4;
sp4[77] = -1 * sp4[76];
sp4[78] = sp4[75] + sp4[77];
sp4[79] = sp4[78] / sp4[16];
sp4[80] = J_c2 * J_c3;
sp4[81] = J_c0 * J_c5;
sp4[82] = -1 * sp4[81];
sp4[83] = sp4[80] + sp4[82];
sp4[84] = sp4[83] / sp4[16];
sp4[85] = J_c0 * J_c4;
sp4[86] = J_c1 * J_c3;
sp4[87] = -1 * sp4[86];
sp4[88] = sp4[85] + sp4[87];
sp4[89] = sp4[88] / sp4[16];
sp4[90] = sp4[84] * sp4[22];
sp4[91] = sp4[89] * sp4[22];
sp4[92] = sp4[79] * sp4[22];
sp4[93] = sp4[84] * sp4[84];
sp4[94] = sp4[84] * sp4[23];
sp4[95] = sp4[84] * sp4[18];
sp4[96] = sp4[89] * sp4[84];
sp4[97] = sp4[79] * sp4[84];
sp4[98] = sp4[89] * sp4[23];
sp4[99] = sp4[79] * sp4[23];
sp4[100] = sp4[89] * sp4[18];
sp4[101] = sp4[89] * sp4[89];
sp4[102] = sp4[89] * sp4[79];
sp4[103] = sp4[79] * sp4[18];
sp4[104] = sp4[79] * sp4[79];
sp4[105] = sp4[69] + sp4[104];
sp4[106] = sp4[70] + sp4[97];
sp4[107] = sp4[71] + sp4[102];
sp4[108] = sp4[72] + sp4[93];
sp4[109] = sp4[73] + sp4[96];
sp4[110] = sp4[74] + sp4[101];
sp4[111] = sp4[37] + sp4[37];
sp4[112] = sp4[42] + sp4[42];
sp4[113] = sp4[47] + sp4[47];
sp4[114] = sp4[112] * sp4[112];
sp4[115] = sp4[112] * sp4[113];
sp4[116] = sp4[111] * sp4[112];
sp4[117] = sp4[113] * sp4[113];
sp4[118] = sp4[111] * sp4[113];
sp4[119] = sp4[111] * sp4[111];
sp4[120] = sp4[119] + sp4[68];
sp4[121] = sp4[116] + sp4[49];
sp4[122] = sp4[118] + sp4[60];
sp4[123] = sp4[114] + sp4[50];
sp4[124] = sp4[115] + sp4[48];
sp4[125] = sp4[117] + sp4[59];
sp4[126] = sp4[42] * sp4[84];
sp4[127] = sp4[84] * sp4[47];
sp4[128] = sp4[37] * sp4[84];
sp4[129] = sp4[42] * sp4[89];
sp4[130] = sp4[42] * sp4[79];
sp4[131] = sp4[89] * sp4[47];
sp4[132] = sp4[37] * sp4[89];
sp4[133] = sp4[79] * sp4[47];
sp4[134] = sp4[37] * sp4[79];
sp4[135] = sp4[120] + sp4[104];
sp4[136] = sp4[121] + sp4[97];
sp4[137] = sp4[122] + sp4[102];
sp4[138] = sp4[123] + sp4[93];
sp4[139] = sp4[124] + sp4[96];
sp4[140] = sp4[125] + sp4[101];
sp4[141] = sp4[105] + sp4[67];
sp4[142] = sp4[106] + sp4[58];
sp4[143] = sp4[107] + sp4[65];
sp4[144] = sp4[66] + sp4[66];
sp4[145] = sp4[53] + sp4[53];
sp4[146] = sp4[62] + sp4[62];
sp4[147] = sp4[108] + sp4[56];
sp4[148] = sp4[109] + sp4[57];
sp4[149] = sp4[55] + sp4[55];
sp4[150] = sp4[51] + sp4[51];
sp4[151] = sp4[54] + sp4[54];
sp4[152] = sp4[110] + sp4[64];
sp4[153] = sp4[63] + sp4[63];
sp4[154] = sp4[52] + sp4[52];
sp4[155] = sp4[61] + sp4[61];
sp4[156] = sp4[135] + sp4[68];
sp4[157] = sp4[136] + sp4[49];
sp4[158] = sp4[137] + sp4[60];
sp4[159] = sp4[138] + sp4[50];
sp4[160] = sp4[139] + sp4[48];
sp4[161] = sp4[140] + sp4[59];
sp4[162] = sp4[67] + sp4[68];
sp4[163] = sp4[58] + sp4[49];
sp4[164] = sp4[65] + sp4[60];
sp4[165] = sp4[56] + sp4[50];
sp4[166] = sp4[57] + sp4[48];
sp4[167] = sp4[59] + sp4[64];
sp4[168] = sp4[79] + sp4[79];
sp4[169] = sp4[84] + sp4[84];
sp4[170] = sp4[89] + sp4[89];
sp4[171] = sp4[168] * sp4[169];
sp4[172] = sp4[169] * sp4[169];
sp4[173] = sp4[170] * sp4[169];
sp4[174] = sp4[170] * sp4[168];
sp4[175] = sp4[170] * sp4[170];
sp4[176] = sp4[168] * sp4[168];
sp4[177] = sp4[162] + sp4[176];
sp4[178] = sp4[163] + sp4[171];
sp4[179] = sp4[164] + sp4[174];
sp4[180] = sp4[165] + sp4[172];
sp4[181] = sp4[166] + sp4[173];
sp4[182] = sp4[167] + sp4[175];
sp4[183] = sp4[141] + sp4[104];
sp4[184] = sp4[142] + sp4[97];
sp4[185] = sp4[143] + sp4[102];
sp4[186] = sp4[103] + sp4[103];
sp4[187] = sp4[92] + sp4[92];
sp4[188] = sp4[99] + sp4[99];
sp4[189] = sp4[147] + sp4[93];
sp4[190] = sp4[148] + sp4[96];
sp4[191] = sp4[95] + sp4[95];
sp4[192] = sp4[90] + sp4[90];
sp4[193] = sp4[94] + sp4[94];
sp4[194] = sp4[152] + sp4[101];
sp4[195] = sp4[100] + sp4[100];
sp4[196] = sp4[91] + sp4[91];
sp4[197] = sp4[98] + sp4[98];
sp4[198] = sp4[156] + sp4[104];
sp4[199] = sp4[157] + sp4[97];
sp4[200] = sp4[158] + sp4[102];
sp4[201] = sp4[134] + sp4[134];
sp4[202] = sp4[130] + sp4[130];
sp4[203] = sp4[133] + sp4[133];
sp4[204] = sp4[159] + sp4[93];
sp4[205] = sp4[160] + sp4[96];
sp4[206] = sp4[128] + sp4[128];
sp4[207] = sp4[126] + sp4[126];
sp4[208] = sp4[127] + sp4[127];
sp4[209] = sp4[161] + sp4[101];
sp4[210] = sp4[132] + sp4[132];
sp4[211] = sp4[129] + sp4[129];
sp4[212] = sp4[131] + sp4[131];
sp4[213] = sp4[177] + sp4[162];
sp4[214] = sp4[178] + sp4[163];
sp4[215] = sp4[179] + sp4[164];
sp4[216] = sp4[180] + sp4[165];
sp4[217] = sp4[181] + sp4[166];
sp4[218] = sp4[182] + sp4[167];
sp4[219] = 0.25 * sp4[191];
sp4[220] = 0.25 * sp4[189];
sp4[221] = 0.25 * sp4[211];
sp4[222] = 0.25 * sp4[216];
sp4[223] = 0.25 * sp4[215];
sp4[224] = 0.25 * sp4[201];
sp4[225] = 0.25 * sp4[187];
sp4[226] = 0.25 * sp4[217];
sp4[227] = 0.25 * sp4[153];
sp4[228] = 0.25 * sp4[154];
sp4[229] = 0.25 * sp4[212];
sp4[230] = 0.25 * sp4[196];
sp4[231] = 0.25 * sp4[200];
sp4[232] = 0.25 * sp4[150];
sp4[233] = 0.25 * sp4[214];
sp4[234] = 0.25 * sp4[199];
sp4[235] = 0.25 * sp4[203];
sp4[236] = 0.25 * sp4[207];
sp4[237] = 0.25 * sp4[146];
sp4[238] = 0.25 * sp4[184];
sp4[239] = 0.25 * sp4[209];
sp4[240] = 0.25 * sp4[197];
sp4[241] = 0.25 * sp4[144];
sp4[242] = 0.25 * sp4[155];
sp4[243] = 0.25 * sp4[188];
sp4[244] = 0.25 * sp4[213];
sp4[245] = 0.25 * sp4[205];
sp4[246] = 0.25 * sp4[193];
sp4[247] = 0.25 * sp4[206];
sp4[248] = 0.25 * sp4[149];
sp4[249] = 0.25 * sp4[202];
sp4[250] = 0.25 * sp4[208];
sp4[251] = 0.25 * sp4[185];
sp4[252] = 0.25 * sp4[192];
sp4[253] = 0.25 * sp4[183];
sp4[254] = 0.25 * sp4[190];
sp4[255] = 0.25 * sp4[145];
sp4[256] = 0.25 * sp4[198];
sp4[257] = 0.25 * sp4[195];
sp4[258] = 0.25 * sp4[218];
sp4[259] = 0.25 * sp4[210];
sp4[260] = 0.25 * sp4[186];
sp4[261] = 0.25 * sp4[204];
sp4[262] = 0.25 * sp4[151];
sp4[263] = 0.25 * sp4[194];
for (int iq = 0; iq < 4; ++iq)
{
    // Quadrature loop body setup (num_points=4)
    // Section for geometrically varying computations
    double w_0 = 0.0;
    for (int ic = 0; ic < 4; ++ic)
    {
        w_0 += w0[0][ic] * FE1_C0_Q4[0][iq][ic - 0];
    }
    double sv4[91];
    sv4[0] = weights4[iq] * sp4[17];
    sv4[1] = sp4[219] * w_0;
    sv4[2] = sp4[220] * w_0;
    sv4[3] = sp4[248] * w_0;
    sv4[4] = sp4[257] * w_0;
    sv4[5] = sp4[224] * w_0;
    sv4[6] = sp4[225] * w_0;
    sv4[7] = sp4[227] * w_0;
    sv4[8] = sp4[228] * w_0;
    sv4[9] = sp4[230] * w_0;
    sv4[10] = sp4[232] * w_0;
    sv4[11] = sp4[234] * w_0;
    sv4[12] = sp4[235] * w_0;
    sv4[13] = sp4[236] * w_0;
    sv4[14] = sp4[238] * w_0;
    sv4[15] = sp4[239] * w_0;
    sv4[16] = sp4[237] * w_0;
    sv4[17] = sp4[243] * w_0;
    sv4[18] = sp4[241] * w_0;
    sv4[19] = sp4[242] * w_0;
    sv4[20] = sp4[221] * w_0;
    sv4[21] = sp4[231] * w_0;
    sv4[22] = sp4[256] * w_0;
    sv4[23] = sp4[244] * w_0;
    sv4[24] = sp4[229] * w_0;
    sv4[25] = sp4[240] * w_0;
    sv4[26] = sp4[258] * w_0;
    sv4[27] = sp4[222] * w_0;
    sv4[28] = sp4[251] * w_0;
    sv4[29] = sp4[226] * w_0;
    sv4[30] = sp4[246] * w_0;
    sv4[31] = sp4[249] * w_0;
    sv4[32] = sp4[250] * w_0;
    sv4[33] = sp4[252] * w_0;
    sv4[34] = sp4[233] * w_0;
    sv4[35] = sp4[255] * w_0;
    sv4[36] = sp4[247] * w_0;
    sv4[37] = sp4[253] * w_0;
    sv4[38] = sp4[223] * w_0;
    sv4[39] = sp4[259] * w_0;
    sv4[40] = sp4[254] * w_0;
    sv4[41] = sp4[260] * w_0;
    sv4[42] = sp4[261] * w_0;
    sv4[43] = sp4[245] * w_0;
    sv4[44] = sp4[262] * w_0;
    sv4[45] = sp4[263] * w_0;
    sv4[46] = sv4[0] * sv4[1];
    sv4[47] = sv4[0] * sv4[2];
    sv4[48] = sv4[0] * sv4[20];
    sv4[49] = sv4[0] * sv4[38];
    sv4[50] = sv4[0] * sv4[5];
    sv4[51] = sv4[0] * sv4[6];
    sv4[52] = sv4[0] * sv4[7];
    sv4[53] = sv4[0] * sv4[8];
    sv4[54] = sv4[0] * sv4[9];
    sv4[55] = sv4[0] * sv4[21];
    sv4[56] = sv4[0] * sv4[10];
    sv4[57] = sv4[0] * sv4[11];
    sv4[58] = sv4[0] * sv4[12];
    sv4[59] = sv4[0] * sv4[13];
    sv4[60] = sv4[0] * sv4[14];
    sv4[61] = sv4[0] * sv4[15];
    sv4[62] = sv4[0] * sv4[16];
    sv4[63] = sv4[0] * sv4[40];
    sv4[64] = sv4[0] * sv4[18];
    sv4[65] = sv4[0] * sv4[19];
    sv4[66] = sv4[0] * sv4[17];
    sv4[67] = sv4[0] * sv4[23];
    sv4[68] = sv4[0] * sv4[24];
    sv4[69] = sv4[0] * sv4[25];
    sv4[70] = sv4[0] * sv4[43];
    sv4[71] = sv4[0] * sv4[27];
    sv4[72] = sv4[0] * sv4[37];
    sv4[73] = sv4[0] * sv4[29];
    sv4[74] = sv4[0] * sv4[30];
    sv4[75] = sv4[0] * sv4[31];
    sv4[76] = sv4[0] * sv4[32];
    sv4[77] = sv4[0] * sv4[28];
    sv4[78] = sv4[0] * sv4[33];
    sv4[79] = sv4[0] * sv4[3];
    sv4[80] = sv4[0] * sv4[34];
    sv4[81] = sv4[0] * sv4[35];
    sv4[82] = sv4[0] * sv4[36];
    sv4[83] = sv4[0] * sv4[22];
    sv4[84] = sv4[0] * sv4[4];
    sv4[85] = sv4[0] * sv4[26];
    sv4[86] = sv4[0] * sv4[39];
    sv4[87] = sv4[0] * sv4[41];
    sv4[88] = sv4[0] * sv4[42];
    sv4[89] = sv4[0] * sv4[44];
    sv4[90] = sv4[0] * sv4[45];
    for (int ia0 = 0; ia0 < 2; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[72] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[60] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[77] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[64] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[81] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[62] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[87] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[51] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[66] * FE0_C0_D100_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 0; ia0 < 3; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[60] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[47] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[63] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[79] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[56] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[89] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[46] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[78] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[74] * FE0_C0_D010_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 0; ia0 < 4; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[77] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[63] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[90] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[52] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[53] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[65] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[84] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[54] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[69] * FE0_C0_D001_Q4[0][iq][ia0 - 0] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 6; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[64] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[79] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[52] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[83] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[57] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[55] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[50] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[75] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[58] * FE0_C0_D100_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 7; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[81] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[56] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[53] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[57] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[88] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[70] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[82] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[59] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[76] * FE0_C0_D010_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 8; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[62] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[89] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[65] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[55] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[70] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[61] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[86] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[48] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[68] * FE0_C0_D001_Q4[0][iq][ia0 - 4] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 10; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[87] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[46] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[84] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[50] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[82] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[86] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[67] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[80] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[49] * FE0_C0_D100_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 11; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[51] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[78] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[54] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[75] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[59] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[48] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[80] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[71] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[73] * FE0_C0_D010_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 12; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv4[66] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv4[74] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv4[69] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv4[58] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv4[76] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv4[68] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv4[49] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D100_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv4[73] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D010_Q4[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv4[85] * FE0_C0_D001_Q4[0][iq][ia0 - 8] * FE0_C0_D001_Q4[0][iq][ia1 - 8];
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
    double buffer_arg0_0[12][12] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec, arg2_0_vec);
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
        
        