
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


void form_cell_integral_0_otherwise(double A[12][12], double *coordinate_dofs, double **w0, double **w1)
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
double sp5[264];
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
sp5[23] = sp5[17] + sp5[17];
sp5[24] = sp5[21] + sp5[21];
sp5[25] = sp5[22] + sp5[22];
sp5[26] = sp5[24] * sp5[24];
sp5[27] = sp5[25] * sp5[24];
sp5[28] = sp5[23] * sp5[24];
sp5[29] = sp5[25] * sp5[25];
sp5[30] = sp5[25] * sp5[23];
sp5[31] = sp5[23] * sp5[23];
sp5[32] = J_c2 * J_c7;
sp5[33] = -1 * J_c1;
sp5[34] = J_c8 * sp5[33];
sp5[35] = sp5[32] + sp5[34];
sp5[36] = sp5[35] / sp5[16];
sp5[37] = J_c0 * J_c8;
sp5[38] = -1 * J_c2;
sp5[39] = J_c6 * sp5[38];
sp5[40] = sp5[37] + sp5[39];
sp5[41] = sp5[40] / sp5[16];
sp5[42] = J_c1 * J_c6;
sp5[43] = J_c0 * J_c7;
sp5[44] = -1 * sp5[43];
sp5[45] = sp5[42] + sp5[44];
sp5[46] = sp5[45] / sp5[16];
sp5[47] = sp5[22] * sp5[22];
sp5[48] = sp5[22] * sp5[46];
sp5[49] = sp5[36] * sp5[22];
sp5[50] = sp5[22] * sp5[17];
sp5[51] = sp5[41] * sp5[22];
sp5[52] = sp5[22] * sp5[21];
sp5[53] = sp5[46] * sp5[17];
sp5[54] = sp5[36] * sp5[17];
sp5[55] = sp5[17] * sp5[17];
sp5[56] = sp5[41] * sp5[17];
sp5[57] = sp5[17] * sp5[21];
sp5[58] = sp5[36] * sp5[46];
sp5[59] = sp5[36] * sp5[36];
sp5[60] = sp5[36] * sp5[41];
sp5[61] = sp5[36] * sp5[21];
sp5[62] = sp5[46] * sp5[21];
sp5[63] = sp5[41] * sp5[21];
sp5[64] = sp5[21] * sp5[21];
sp5[65] = sp5[41] * sp5[46];
sp5[66] = sp5[41] * sp5[41];
sp5[67] = sp5[46] * sp5[46];
sp5[68] = sp5[31] + sp5[59];
sp5[69] = sp5[28] + sp5[60];
sp5[70] = sp5[30] + sp5[58];
sp5[71] = sp5[26] + sp5[66];
sp5[72] = sp5[27] + sp5[65];
sp5[73] = sp5[29] + sp5[67];
sp5[74] = J_c1 * J_c5;
sp5[75] = J_c2 * J_c4;
sp5[76] = -1 * sp5[75];
sp5[77] = sp5[74] + sp5[76];
sp5[78] = sp5[77] / sp5[16];
sp5[79] = J_c2 * J_c3;
sp5[80] = J_c0 * J_c5;
sp5[81] = -1 * sp5[80];
sp5[82] = sp5[79] + sp5[81];
sp5[83] = sp5[82] / sp5[16];
sp5[84] = J_c0 * J_c4;
sp5[85] = J_c1 * J_c3;
sp5[86] = -1 * sp5[85];
sp5[87] = sp5[84] + sp5[86];
sp5[88] = sp5[87] / sp5[16];
sp5[89] = sp5[88] * sp5[22];
sp5[90] = sp5[78] * sp5[22];
sp5[91] = sp5[83] * sp5[22];
sp5[92] = sp5[88] * sp5[17];
sp5[93] = sp5[78] * sp5[17];
sp5[94] = sp5[83] * sp5[17];
sp5[95] = sp5[88] * sp5[78];
sp5[96] = sp5[78] * sp5[78];
sp5[97] = sp5[78] * sp5[21];
sp5[98] = sp5[78] * sp5[83];
sp5[99] = sp5[88] * sp5[83];
sp5[100] = sp5[83] * sp5[21];
sp5[101] = sp5[83] * sp5[83];
sp5[102] = sp5[88] * sp5[88];
sp5[103] = sp5[88] * sp5[21];
sp5[104] = sp5[68] + sp5[96];
sp5[105] = sp5[69] + sp5[98];
sp5[106] = sp5[70] + sp5[95];
sp5[107] = sp5[71] + sp5[101];
sp5[108] = sp5[72] + sp5[99];
sp5[109] = sp5[73] + sp5[102];
sp5[110] = sp5[36] + sp5[36];
sp5[111] = sp5[41] + sp5[41];
sp5[112] = sp5[46] + sp5[46];
sp5[113] = sp5[111] * sp5[112];
sp5[114] = sp5[110] * sp5[111];
sp5[115] = sp5[111] * sp5[111];
sp5[116] = sp5[112] * sp5[112];
sp5[117] = sp5[110] * sp5[112];
sp5[118] = sp5[110] * sp5[110];
sp5[119] = sp5[118] + sp5[55];
sp5[120] = sp5[114] + sp5[57];
sp5[121] = sp5[117] + sp5[50];
sp5[122] = sp5[115] + sp5[64];
sp5[123] = sp5[113] + sp5[52];
sp5[124] = sp5[116] + sp5[47];
sp5[125] = sp5[88] * sp5[46];
sp5[126] = sp5[36] * sp5[88];
sp5[127] = sp5[41] * sp5[88];
sp5[128] = sp5[78] * sp5[46];
sp5[129] = sp5[36] * sp5[78];
sp5[130] = sp5[41] * sp5[78];
sp5[131] = sp5[83] * sp5[46];
sp5[132] = sp5[36] * sp5[83];
sp5[133] = sp5[41] * sp5[83];
sp5[134] = sp5[119] + sp5[96];
sp5[135] = sp5[120] + sp5[98];
sp5[136] = sp5[121] + sp5[95];
sp5[137] = sp5[122] + sp5[101];
sp5[138] = sp5[123] + sp5[99];
sp5[139] = sp5[124] + sp5[102];
sp5[140] = sp5[104] + sp5[59];
sp5[141] = sp5[105] + sp5[60];
sp5[142] = sp5[106] + sp5[58];
sp5[143] = sp5[54] + sp5[54];
sp5[144] = sp5[61] + sp5[61];
sp5[145] = sp5[49] + sp5[49];
sp5[146] = sp5[107] + sp5[66];
sp5[147] = sp5[108] + sp5[65];
sp5[148] = sp5[56] + sp5[56];
sp5[149] = sp5[63] + sp5[63];
sp5[150] = sp5[51] + sp5[51];
sp5[151] = sp5[109] + sp5[67];
sp5[152] = sp5[53] + sp5[53];
sp5[153] = sp5[62] + sp5[62];
sp5[154] = sp5[48] + sp5[48];
sp5[155] = sp5[134] + sp5[55];
sp5[156] = sp5[135] + sp5[57];
sp5[157] = sp5[136] + sp5[50];
sp5[158] = sp5[137] + sp5[64];
sp5[159] = sp5[138] + sp5[52];
sp5[160] = sp5[139] + sp5[47];
sp5[161] = sp5[59] + sp5[55];
sp5[162] = sp5[60] + sp5[57];
sp5[163] = sp5[58] + sp5[50];
sp5[164] = sp5[66] + sp5[64];
sp5[165] = sp5[65] + sp5[52];
sp5[166] = sp5[47] + sp5[67];
sp5[167] = sp5[78] + sp5[78];
sp5[168] = sp5[83] + sp5[83];
sp5[169] = sp5[88] + sp5[88];
sp5[170] = sp5[168] * sp5[168];
sp5[171] = sp5[169] * sp5[168];
sp5[172] = sp5[167] * sp5[168];
sp5[173] = sp5[169] * sp5[169];
sp5[174] = sp5[169] * sp5[167];
sp5[175] = sp5[167] * sp5[167];
sp5[176] = sp5[161] + sp5[175];
sp5[177] = sp5[162] + sp5[172];
sp5[178] = sp5[163] + sp5[174];
sp5[179] = sp5[164] + sp5[170];
sp5[180] = sp5[165] + sp5[171];
sp5[181] = sp5[166] + sp5[173];
sp5[182] = sp5[140] + sp5[96];
sp5[183] = sp5[141] + sp5[98];
sp5[184] = sp5[142] + sp5[95];
sp5[185] = sp5[93] + sp5[93];
sp5[186] = sp5[97] + sp5[97];
sp5[187] = sp5[90] + sp5[90];
sp5[188] = sp5[146] + sp5[101];
sp5[189] = sp5[147] + sp5[99];
sp5[190] = sp5[94] + sp5[94];
sp5[191] = sp5[100] + sp5[100];
sp5[192] = sp5[91] + sp5[91];
sp5[193] = sp5[151] + sp5[102];
sp5[194] = sp5[92] + sp5[92];
sp5[195] = sp5[103] + sp5[103];
sp5[196] = sp5[89] + sp5[89];
sp5[197] = sp5[155] + sp5[96];
sp5[198] = sp5[156] + sp5[98];
sp5[199] = sp5[157] + sp5[95];
sp5[200] = sp5[129] + sp5[129];
sp5[201] = sp5[130] + sp5[130];
sp5[202] = sp5[128] + sp5[128];
sp5[203] = sp5[158] + sp5[101];
sp5[204] = sp5[159] + sp5[99];
sp5[205] = sp5[132] + sp5[132];
sp5[206] = sp5[133] + sp5[133];
sp5[207] = sp5[131] + sp5[131];
sp5[208] = sp5[160] + sp5[102];
sp5[209] = sp5[126] + sp5[126];
sp5[210] = sp5[127] + sp5[127];
sp5[211] = sp5[125] + sp5[125];
sp5[212] = sp5[176] + sp5[161];
sp5[213] = sp5[177] + sp5[162];
sp5[214] = sp5[178] + sp5[163];
sp5[215] = sp5[179] + sp5[164];
sp5[216] = sp5[180] + sp5[165];
sp5[217] = sp5[181] + sp5[166];
sp5[218] = 0.25 * sp5[216];
sp5[219] = 0.25 * sp5[211];
sp5[220] = 0.25 * sp5[191];
sp5[221] = 0.25 * sp5[199];
sp5[222] = 0.25 * sp5[204];
sp5[223] = 0.25 * sp5[145];
sp5[224] = 0.25 * sp5[184];
sp5[225] = 0.25 * sp5[150];
sp5[226] = 0.25 * sp5[198];
sp5[227] = 0.25 * sp5[153];
sp5[228] = 0.25 * sp5[207];
sp5[229] = 0.25 * sp5[144];
sp5[230] = 0.25 * sp5[203];
sp5[231] = 0.25 * sp5[186];
sp5[232] = 0.25 * sp5[154];
sp5[233] = 0.25 * sp5[195];
sp5[234] = 0.25 * sp5[206];
sp5[235] = 0.25 * sp5[148];
sp5[236] = 0.25 * sp5[205];
sp5[237] = 0.25 * sp5[183];
sp5[238] = 0.25 * sp5[196];
sp5[239] = 0.25 * sp5[152];
sp5[240] = 0.25 * sp5[209];
sp5[241] = 0.25 * sp5[214];
sp5[242] = 0.25 * sp5[192];
sp5[243] = 0.25 * sp5[200];
sp5[244] = 0.25 * sp5[143];
sp5[245] = 0.25 * sp5[194];
sp5[246] = 0.25 * sp5[201];
sp5[247] = 0.25 * sp5[189];
sp5[248] = 0.25 * sp5[149];
sp5[249] = 0.25 * sp5[210];
sp5[250] = 0.25 * sp5[197];
sp5[251] = 0.25 * sp5[193];
sp5[252] = 0.25 * sp5[190];
sp5[253] = 0.25 * sp5[202];
sp5[254] = 0.25 * sp5[213];
sp5[255] = 0.25 * sp5[185];
sp5[256] = 0.25 * sp5[187];
sp5[257] = 0.25 * sp5[217];
sp5[258] = 0.25 * sp5[215];
sp5[259] = 0.25 * sp5[208];
sp5[260] = 0.25 * sp5[212];
sp5[261] = 0.25 * sp5[188];
sp5[262] = 0.25 * sp5[182];
sp5[263] = fabs(sp5[16]);
for (int iq = 0; iq < 5; ++iq)
{
    // Quadrature loop body setup (num_points=5)
    // Section for geometrically varying computations
    double w_0 = 0.0;
    for (int ic = 0; ic < 4; ++ic)
    {
        w_0 += w0[0][ic] * FE1_C0_Q5[0][iq][ic - 0];
    }
    double w_1 = 0.0;
    for (int ic = 0; ic < 4; ++ic)
    {
        w_1 += w1[0][ic] * FE1_C0_Q5[0][iq][ic - 0];
    }
    double sv5[92];
    sv5[0] = w_0 * w_1;
    sv5[1] = sp5[218] * sv5[0];
    sv5[2] = sp5[220] * sv5[0];
    sv5[3] = sp5[257] * sv5[0];
    sv5[4] = sp5[222] * sv5[0];
    sv5[5] = sp5[223] * sv5[0];
    sv5[6] = sp5[224] * sv5[0];
    sv5[7] = sp5[246] * sv5[0];
    sv5[8] = sp5[225] * sv5[0];
    sv5[9] = sp5[226] * sv5[0];
    sv5[10] = sp5[227] * sv5[0];
    sv5[11] = sp5[228] * sv5[0];
    sv5[12] = sp5[229] * sv5[0];
    sv5[13] = sp5[231] * sv5[0];
    sv5[14] = sp5[232] * sv5[0];
    sv5[15] = sp5[219] * sv5[0];
    sv5[16] = sp5[233] * sv5[0];
    sv5[17] = sp5[234] * sv5[0];
    sv5[18] = sp5[235] * sv5[0];
    sv5[19] = sp5[236] * sv5[0];
    sv5[20] = sp5[237] * sv5[0];
    sv5[21] = sp5[238] * sv5[0];
    sv5[22] = sp5[239] * sv5[0];
    sv5[23] = sp5[240] * sv5[0];
    sv5[24] = sp5[221] * sv5[0];
    sv5[25] = sp5[243] * sv5[0];
    sv5[26] = sp5[242] * sv5[0];
    sv5[27] = sp5[230] * sv5[0];
    sv5[28] = sp5[245] * sv5[0];
    sv5[29] = sp5[256] * sv5[0];
    sv5[30] = sp5[247] * sv5[0];
    sv5[31] = sp5[241] * sv5[0];
    sv5[32] = sp5[249] * sv5[0];
    sv5[33] = sp5[250] * sv5[0];
    sv5[34] = sp5[251] * sv5[0];
    sv5[35] = sp5[252] * sv5[0];
    sv5[36] = sp5[253] * sv5[0];
    sv5[37] = sp5[255] * sv5[0];
    sv5[38] = sp5[244] * sv5[0];
    sv5[39] = sp5[248] * sv5[0];
    sv5[40] = sp5[258] * sv5[0];
    sv5[41] = sp5[259] * sv5[0];
    sv5[42] = sp5[260] * sv5[0];
    sv5[43] = sp5[254] * sv5[0];
    sv5[44] = sp5[261] * sv5[0];
    sv5[45] = sp5[262] * sv5[0];
    sv5[46] = weights5[iq] * sp5[263];
    sv5[47] = sv5[1] * sv5[46];
    sv5[48] = sv5[15] * sv5[46];
    sv5[49] = sv5[2] * sv5[46];
    sv5[50] = sv5[24] * sv5[46];
    sv5[51] = sv5[4] * sv5[46];
    sv5[52] = sv5[5] * sv5[46];
    sv5[53] = sv5[6] * sv5[46];
    sv5[54] = sv5[8] * sv5[46];
    sv5[55] = sv5[9] * sv5[46];
    sv5[56] = sv5[10] * sv5[46];
    sv5[57] = sv5[11] * sv5[46];
    sv5[58] = sv5[12] * sv5[46];
    sv5[59] = sv5[13] * sv5[46];
    sv5[60] = sv5[14] * sv5[46];
    sv5[61] = sv5[16] * sv5[46];
    sv5[62] = sv5[18] * sv5[46];
    sv5[63] = sv5[19] * sv5[46];
    sv5[64] = sv5[20] * sv5[46];
    sv5[65] = sv5[21] * sv5[46];
    sv5[66] = sv5[22] * sv5[46];
    sv5[67] = sv5[25] * sv5[46];
    sv5[68] = sv5[23] * sv5[46];
    sv5[69] = sv5[31] * sv5[46];
    sv5[70] = sv5[7] * sv5[46];
    sv5[71] = sv5[17] * sv5[46];
    sv5[72] = sv5[26] * sv5[46];
    sv5[73] = sv5[27] * sv5[46];
    sv5[74] = sv5[28] * sv5[46];
    sv5[75] = sv5[30] * sv5[46];
    sv5[76] = sv5[3] * sv5[46];
    sv5[77] = sv5[33] * sv5[46];
    sv5[78] = sv5[34] * sv5[46];
    sv5[79] = sv5[35] * sv5[46];
    sv5[80] = sv5[36] * sv5[46];
    sv5[81] = sv5[43] * sv5[46];
    sv5[82] = sv5[37] * sv5[46];
    sv5[83] = sv5[38] * sv5[46];
    sv5[84] = sv5[39] * sv5[46];
    sv5[85] = sv5[29] * sv5[46];
    sv5[86] = sv5[40] * sv5[46];
    sv5[87] = sv5[32] * sv5[46];
    sv5[88] = sv5[41] * sv5[46];
    sv5[89] = sv5[42] * sv5[46];
    sv5[90] = sv5[44] * sv5[46];
    sv5[91] = sv5[45] * sv5[46];
    for (int ia0 = 0; ia0 < 2; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[91] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[64] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[53] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[83] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[58] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[52] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[82] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[59] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[85] * FE0_C0_D100_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 0; ia0 < 3; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[64] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[90] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[75] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[62] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[84] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[54] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[79] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[49] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[72] * FE0_C0_D010_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 0; ia0 < 4; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[53] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[75] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[78] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[66] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[56] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[60] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[74] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[61] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[65] * FE0_C0_D001_Q5[0][iq][ia0 - 0] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 6; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[83] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[62] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[66] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[77] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[55] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[50] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[67] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[70] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[80] * FE0_C0_D100_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 7; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[58] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[84] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[56] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[55] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[73] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[51] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[63] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[71] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[57] * FE0_C0_D010_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 4; ia0 < 8; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[52] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[54] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[60] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[50] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[51] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[88] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[68] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[87] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[48] * FE0_C0_D001_Q5[0][iq][ia0 - 4] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 10; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[82] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[79] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[74] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[67] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[63] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[68] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[89] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[81] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[69] * FE0_C0_D100_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 11; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[59] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[49] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[61] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[70] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[71] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[87] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[81] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[86] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[47] * FE0_C0_D010_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
    for (int ia0 = 8; ia0 < 12; ++ia0)
    {
        for (int ia1 = 0; ia1 < 2; ++ia1)
        {
            A[ia0][ia1] += sv5[85] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 3; ++ia1)
        {
            A[ia0][ia1] += sv5[72] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 0; ia1 < 4; ++ia1)
        {
            A[ia0][ia1] += sv5[65] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 0];
        }
        for (int ia1 = 4; ia1 < 6; ++ia1)
        {
            A[ia0][ia1] += sv5[80] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 7; ++ia1)
        {
            A[ia0][ia1] += sv5[57] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 4; ia1 < 8; ++ia1)
        {
            A[ia0][ia1] += sv5[48] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 4];
        }
        for (int ia1 = 8; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv5[69] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D100_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 11; ++ia1)
        {
            A[ia0][ia1] += sv5[47] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D010_Q5[0][iq][ia1 - 8];
        }
        for (int ia1 = 8; ia1 < 12; ++ia1)
        {
            A[ia0][ia1] += sv5[76] * FE0_C0_D001_Q5[0][iq][ia0 - 8] * FE0_C0_D001_Q5[0][iq][ia1 - 8];
        }
    }
}
}

            

        
        void wrap_form_cell_integral_0_otherwise(int start, int end,
                      Mat arg0_0_, int *arg0_0_map0_0, int *arg0_0_map1_0, double *arg1_0, int *arg1_0_map0_0, double *arg2_0, int *arg2_0_map0_0, double *arg3_0, int *arg3_0_map0_0
                      ) {
  Mat arg0_0_0 = arg0_0_;
  double *arg1_0_vec[12];
    double *arg2_0_vec[4];
    double *arg3_0_vec[4];
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
    arg3_0_vec[0] = arg3_0 + (arg3_0_map0_0[i * 4 + 0])* 1;
    arg3_0_vec[1] = arg3_0 + (arg3_0_map0_0[i * 4 + 1])* 1;
    arg3_0_vec[2] = arg3_0 + (arg3_0_map0_0[i * 4 + 2])* 1;
    arg3_0_vec[3] = arg3_0 + (arg3_0_map0_0[i * 4 + 3])* 1;
    double buffer_arg0_0[12][12] __attribute__((aligned(32))) = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec, arg2_0_vec, arg3_0_vec);
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
        
        