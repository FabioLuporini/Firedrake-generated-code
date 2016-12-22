
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


void form_cell_integral_0_otherwise(double A[10][10], double *coordinate_dofs)
{
// Section for quadrature weights and points
static const double weights15[15] = { 0.0302836780971, 0.00602678571429, 0.00602678571429, 0.00602678571429, 0.00602678571429, 0.011645249086, 0.011645249086, 0.011645249086, 0.011645249086, 0.0109491415614, 0.0109491415614, 0.0109491415614, 0.0109491415614, 0.0109491415614, 0.0109491415614 };
static const double points15[45] = { 0.25, 0.25, 0.25, 0.0, 0.333333333333, 0.333333333333, 0.333333333333, 0.333333333333, 0.333333333333, 0.333333333333, 0.333333333333, 0.0, 0.333333333333, 0.0, 0.333333333333, 0.727272727273, 0.0909090909091, 0.0909090909091, 0.0909090909091, 0.0909090909091, 0.0909090909091, 0.0909090909091, 0.0909090909091, 0.727272727273, 0.0909090909091, 0.727272727273, 0.0909090909091, 0.433449846426, 0.0665501535737, 0.0665501535737, 0.0665501535737, 0.433449846426, 0.0665501535737, 0.0665501535737, 0.0665501535737, 0.433449846426, 0.0665501535737, 0.433449846426, 0.433449846426, 0.433449846426, 0.0665501535737, 0.433449846426, 0.433449846426, 0.433449846426, 0.0665501535737 };
// Section for precomputed element basis function values
// Table dimensions: num_entities, num_points, num_dofs
// Definitions of 7 tables for 15 quadrature points
static const double FE0_C0_D001_Q15[1][15][4] =
    { { { -1.0, -5.55111512313e-17, 0.0, 1.0 },
        { -1.0, -5.55111512313e-17, 0.0, 1.0 },
        { -1.0, -5.55111512313e-17, 0.0, 1.0 },
        { -1.0, 1.11022302463e-16, 5.55111512313e-17, 1.0 },
        { -1.0, -5.55111512313e-17, 0.0, 1.0 },
        { -1.0, 1.11022302463e-16, 5.55111512313e-17, 1.0 },
        { -1.0, 1.11022302463e-16, 5.55111512313e-17, 1.0 },
        { -1.0, -2.77555756156e-16, -1.11022302463e-16, 1.0 },
        { -1.0, 1.66533453694e-16, 1.66533453694e-16, 1.0 },
        { -1.0, 1.11022302463e-16, 5.55111512313e-17, 1.0 },
        { -1.0, 1.11022302463e-16, 5.55111512313e-17, 1.0 },
        { -1.0, -1.66533453694e-16, -5.55111512313e-17, 1.0 },
        { -1.0, -1.66533453694e-16, -5.55111512313e-17, 1.0 },
        { -1.0, -1.66533453694e-16, -5.55111512313e-17, 1.0 },
        { -1.0, 1.11022302463e-16, 5.55111512313e-17, 1.0 } } };
static const double FE0_C0_D010_Q15[1][15][3] =
    { { { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -5.55111512313e-17, 1.0 },
        { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -1.11022302463e-16, 1.0 },
        { -1.0, -2.77555756156e-16, 1.0 },
        { -1.0, 0.0, 1.0 },
        { -1.0, -5.55111512313e-17, 1.0 },
        { -1.0, -5.55111512313e-17, 1.0 },
        { -1.0, -1.66533453694e-16, 1.0 },
        { -1.0, -1.66533453694e-16, 1.0 },
        { -1.0, -1.66533453694e-16, 1.0 },
        { -1.0, -5.55111512313e-17, 1.0 } } };
static const double FE0_C0_D100_Q15[1][15][2] =
    { { { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 },
        { -1.0, 1.0 } } };
static const double FE1_C0_D001_Q15[1][15][10] =
    { { { 4.64905891562e-16, -1.7694179455e-15, 5.55111512313e-16, -3.88578058619e-15, 1.0, 1.0, -2.33146835171e-15, -2.49800180541e-16, -1.0, -1.0 },
        { -0.333333333333, -2.04377742368e-15, 1.06484704527e-15, 0.333333333333, 1.33333333333, 7.56448703119e-15, -2.65440008723e-15, -7.62185653501e-16, -1.33333333333, -4.61709613592e-16 },
        { 1.0, -2.00719087701e-15, 9.6265571386e-16, 0.333333333333, 1.33333333333, 1.33333333333, -3.17292505287e-15, -1.33333333333, -1.33333333333, -1.33333333333 },
        { -0.333333333333, -1.11022302463e-15, -5.55111512313e-17, -1.0, 1.33333333333, 1.33333333333, -3.1918911958e-15, 1.33333333333, -1.33333333333, -1.33333333333 },
        { -0.333333333333, -2.08101992064e-15, 3.61470733537e-16, 0.333333333333, -5.57704514446e-16, 1.33333333333, -7.20996715473e-16, 7.19051963873e-16, 3.61470733537e-16, -1.33333333333 },
        { 0.636363636364, -3.16413562018e-15, -3.05311331772e-16, -0.636363636364, 0.363636363636, 2.90909090909, -2.41473507856e-15, -2.99760216649e-15, -0.363636363636, -2.90909090909 },
        { -1.90909090909, -1.27675647832e-15, -1.38777878078e-16, -0.636363636364, 0.363636363636, 0.363636363636, -1.60982338571e-15, 2.54545454545, -0.363636363636, -0.363636363636 },
        { 0.636363636364, -2.49800180541e-15, 1.38777878078e-15, 1.90909090909, 0.363636363636, 0.363636363636, -8.32667268469e-16, -2.54545454545, -0.363636363636, -0.363636363636 },
        { 0.636363636364, 2.22044604925e-16, 1.27675647832e-15, -0.636363636364, 2.90909090909, 0.363636363636, -5.02375918643e-15, 4.4408920985e-16, -2.90909090909, -0.363636363636 },
        { -0.733799385705, -1.36002320517e-15, -1.94289029309e-16, -0.733799385705, 0.266200614295, 1.73379938571, -1.55431223448e-15, 1.46759877141, -0.266200614295, -1.73379938571 },
        { -0.733799385705, -1.49880108324e-15, 3.88578058619e-16, -0.733799385705, 1.73379938571, 0.266200614295, -3.13638004457e-15, 1.46759877141, -1.73379938571, -0.266200614295 },
        { -0.733799385705, -2.24820162487e-15, 8.18789480661e-16, 0.733799385705, 0.266200614295, 0.266200614295, -1.09634523682e-15, -1.11022302463e-16, -0.266200614295, -0.266200614295 },
        { 0.733799385705, -1.65145674913e-15, 1.29063426613e-15, 0.733799385705, 1.73379938571, 0.266200614295, -3.06699110553e-15, -1.46759877141, -1.73379938571, -0.266200614295 },
        { 0.733799385705, -2.92821322745e-15, 6.80011602583e-16, 0.733799385705, 0.266200614295, 1.73379938571, -1.16573417586e-15, -1.46759877141, -0.266200614295, -1.73379938571 },
        { 0.733799385705, -1.36002320517e-15, 2.49800180541e-16, -0.733799385705, 1.73379938571, 1.73379938571, -4.05231403988e-15, -1.22124532709e-15, -1.73379938571, -1.73379938571 } } };
static const double FE1_C0_D010_Q15[1][15][10] =
    { { { 9.85322934355e-16, -2.53963516883e-15, -2.99760216649e-15, -2.77555756156e-17, 1.0, -3.05311331772e-16, 1.0, -1.0, -6.52256026967e-16, -1.0 },
        { -0.333333333333, -2.47024622979e-15, 0.333333333333, -3.70074341542e-17, 1.33333333333, -1.11022302463e-16, 5.96744875736e-15, -1.33333333333, -1.66533453694e-16, -9.15933995316e-16 },
        { 1.0, -2.80331313718e-15, 0.333333333333, -3.70074341542e-17, 1.33333333333, -2.77555756156e-16, 1.33333333333, -1.33333333333, -1.33333333333, -1.33333333333 },
        { -0.333333333333, -3.33066907388e-15, 0.333333333333, -8.01186856865e-32, 8.881784197e-16, -4.4408920985e-16, 1.33333333333, -4.4408920985e-16, -1.01307850997e-15, -1.33333333333 },
        { -0.333333333333, -1.80411241502e-15, -1.0, -3.70074341542e-17, 1.33333333333, -1.66533453694e-16, 1.33333333333, -1.33333333333, 1.33333333333, -1.33333333333 },
        { 0.636363636364, -3.51801920928e-15, -0.636363636364, -1.00929365875e-17, 0.363636363636, -8.32667268469e-17, 2.90909090909, -0.363636363636, -2.48065457065e-15, -2.90909090909 },
        { -1.90909090909, -2.24126273096e-15, -0.636363636364, -1.00929365875e-17, 0.363636363636, -2.49800180541e-16, 0.363636363636, -0.363636363636, 2.54545454545, -0.363636363636 },
        { 0.636363636364, -1.33226762955e-15, -0.636363636364, -8.07434927e-17, 2.90909090909, 0.0, 0.363636363636, -2.90909090909, -1.66533453694e-15, -0.363636363636 },
        { 0.636363636364, -2.92821322745e-15, 1.90909090909, -1.00929365875e-17, 0.363636363636, -4.4408920985e-16, 0.363636363636, -0.363636363636, -2.54545454545, -0.363636363636 },
        { -0.733799385705, -2.81719092499e-15, -0.733799385705, -7.38855127898e-18, 0.266200614295, -1.94289029309e-16, 1.73379938571, -0.266200614295, 1.46759877141, -1.73379938571 },
        { -0.733799385705, -3.20576898361e-15, 0.733799385705, -7.38855127898e-18, 0.266200614295, -2.22044604925e-16, 0.266200614295, -0.266200614295, -9.2287288922e-16, -0.266200614295 },
        { -0.733799385705, -1.47104550763e-15, -0.733799385705, -4.81225999523e-17, 1.73379938571, -3.33066907388e-16, 0.266200614295, -1.73379938571, 1.46759877141, -0.266200614295 },
        { 0.733799385705, -2.77555756156e-15, 0.733799385705, -4.81225999523e-17, 1.73379938571, -2.22044604925e-16, 0.266200614295, -1.73379938571, -1.46759877141, -0.266200614295 },
        { 0.733799385705, -2.41473507856e-15, -0.733799385705, -4.81225999523e-17, 1.73379938571, -2.22044604925e-16, 1.73379938571, -1.73379938571, -1.72084568817e-15, -1.73379938571 },
        { 0.733799385705, -3.56659146661e-15, 0.733799385705, -7.38855127898e-18, 0.266200614295, -3.60822483003e-16, 1.73379938571, -0.266200614295, -1.46759877141, -1.73379938571 } } };
static const double FE1_C0_D100_Q15[1][15][10] =
    { { { 6.66133814775e-16, -7.77156117238e-16, 4.99600361081e-17, 1.11022302463e-16, -7.39557098645e-33, 1.0, 1.0, -1.0, -1.0, 1.88962527671e-16 },
        { -0.333333333333, -1.0, 6.66133814775e-17, 1.48029736617e-16, -9.86076131526e-33, 1.33333333333, 1.33333333333, -1.33333333333, -1.33333333333, 1.33333333333 },
        { 1.0, 0.333333333333, 6.66133814775e-17, 1.48029736617e-16, -9.86076131526e-33, 1.33333333333, 1.33333333333, -1.33333333333, -1.33333333333, -1.33333333333 },
        { -0.333333333333, 0.333333333333, 1.48029736617e-16, 0.0, -8.21730109605e-33, 0.0, 1.33333333333, 0.0, -1.33333333333, -4.54024556813e-16 },
        { -0.333333333333, 0.333333333333, -8.14163551392e-17, 1.48029736617e-16, -1.64346021921e-33, 1.33333333333, -5.55111512313e-17, -1.33333333333, 0.0, 1.58865213028e-15 },
        { 0.636363636364, 1.90909090909, 1.81672858575e-17, 4.037174635e-17, -2.68929854053e-33, 0.363636363636, 0.363636363636, -0.363636363636, -0.363636363636, -2.54545454545 },
        { -1.90909090909, -0.636363636364, 1.81672858575e-17, 4.037174635e-17, -2.68929854053e-33, 0.363636363636, 0.363636363636, -0.363636363636, -0.363636363636, 2.54545454545 },
        { 0.636363636364, -0.636363636364, -1.3726393759e-16, 3.229739708e-16, -5.82681350447e-33, 2.90909090909, 0.363636363636, -2.90909090909, -0.363636363636, -1.79864961371e-15 },
        { 0.636363636364, -0.636363636364, 3.00769510308e-16, 4.037174635e-17, -1.83768733603e-32, 0.363636363636, 2.90909090909, -0.363636363636, -2.90909090909, 2.01650877335e-15 },
        { -0.733799385705, 0.733799385705, 1.32993923022e-17, 2.95542051159e-17, -1.96870553965e-33, 0.266200614295, 0.266200614295, -0.266200614295, -0.266200614295, -4.2446866458e-16 },
        { -0.733799385705, -0.733799385705, 1.76235586995e-16, 2.95542051159e-17, -1.10134812843e-32, 0.266200614295, 1.73379938571, -0.266200614295, -1.73379938571, 1.46759877141 },
        { -0.733799385705, -0.733799385705, -7.63155147791e-17, 1.92490399809e-16, -3.77766068858e-33, 1.73379938571, 0.266200614295, -1.73379938571, -0.266200614295, 1.46759877141 },
        { 0.733799385705, -0.733799385705, 8.66206799141e-17, 1.92490399809e-16, -1.28224364332e-32, 1.73379938571, 1.73379938571, -1.73379938571, -1.73379938571, -1.11556759811e-15 },
        { 0.733799385705, 0.733799385705, -7.63155147791e-17, 1.92490399809e-16, -3.77766068858e-33, 1.73379938571, 0.266200614295, -1.73379938571, -0.266200614295, -1.46759877141 },
        { 0.733799385705, 0.733799385705, 1.76235586995e-16, 2.95542051159e-17, -1.10134812843e-32, 0.266200614295, 1.73379938571, -0.266200614295, -1.73379938571, -1.46759877141 } } };
static const double FE1_C0_Q15[1][15][10] =
    { { { -0.125, -0.125, -0.125, -0.125, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 },
        { -0.111111111111, 1.21430643318e-17, -0.111111111111, -0.111111111111, 0.444444444444, 4.16333634234e-17, 4.85722573274e-17, 0.444444444444, 0.444444444444, -2.60208521397e-17 },
        { -6.07153216592e-17, -0.111111111111, -0.111111111111, -0.111111111111, 0.444444444444, 0.444444444444, 0.444444444444, 8.32667268469e-17, -1.17961196366e-16, 5.72458747072e-17 },
        { -0.111111111111, -0.111111111111, -0.111111111111, -2.77555756156e-17, 1.38777878078e-17, 1.38777878078e-17, 0.444444444444, 5.55111512313e-17, 0.444444444444, 0.444444444444 },
        { -0.111111111111, -0.111111111111, 0.0, -0.111111111111, -8.32667268469e-17, 0.444444444444, -2.77555756156e-17, 0.444444444444, -5.37764277553e-17, 0.444444444444 },
        { -0.0743801652893, 0.330578512397, -0.0743801652893, -0.0743801652893, 0.0330578512397, 0.264462809917, 0.264462809917, 0.0330578512397, 0.0330578512397, 0.264462809917 },
        { 0.330578512397, -0.0743801652893, -0.0743801652893, -0.0743801652893, 0.0330578512397, 0.0330578512397, 0.0330578512397, 0.264462809917, 0.264462809917, 0.264462809917 },
        { -0.0743801652893, -0.0743801652893, -0.0743801652893, 0.330578512397, 0.264462809917, 0.264462809917, 0.0330578512397, 0.264462809917, 0.0330578512397, 0.0330578512397 },
        { -0.0743801652893, -0.0743801652893, 0.330578512397, -0.0743801652893, 0.264462809917, 0.0330578512397, 0.264462809917, 0.0330578512397, 0.264462809917, 0.0330578512397 },
        { -0.0576923076923, -0.0576923076923, -0.0576923076923, -0.0576923076923, 0.0177156917627, 0.115384615385, 0.115384615385, 0.115384615385, 0.115384615385, 0.751515077468 },
        { -0.0576923076923, -0.0576923076923, -0.0576923076923, -0.0576923076923, 0.115384615385, 0.0177156917627, 0.115384615385, 0.115384615385, 0.751515077468, 0.115384615385 },
        { -0.0576923076923, -0.0576923076923, -0.0576923076923, -0.0576923076923, 0.115384615385, 0.115384615385, 0.0177156917627, 0.751515077468, 0.115384615385, 0.115384615385 },
        { -0.0576923076923, -0.0576923076923, -0.0576923076923, -0.0576923076923, 0.751515077468, 0.115384615385, 0.115384615385, 0.115384615385, 0.115384615385, 0.0177156917627 },
        { -0.0576923076923, -0.0576923076923, -0.0576923076923, -0.0576923076923, 0.115384615385, 0.751515077468, 0.115384615385, 0.115384615385, 0.0177156917627, 0.115384615385 },
        { -0.0576923076923, -0.0576923076923, -0.0576923076923, -0.0576923076923, 0.115384615385, 0.115384615385, 0.751515077468, 0.0177156917627, 0.115384615385, 0.115384615385 } } };
// Section for piecewise constant computations
const double J_c4 = coordinate_dofs[1] * FE0_C0_D010_Q15[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D010_Q15[0][0][5 - 4] + coordinate_dofs[7] * FE0_C0_D010_Q15[0][0][6 - 4];
const double J_c8 = coordinate_dofs[2] * FE0_C0_D001_Q15[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D001_Q15[0][0][9 - 8] + coordinate_dofs[8] * FE0_C0_D001_Q15[0][0][10 - 8] + coordinate_dofs[11] * FE0_C0_D001_Q15[0][0][11 - 8];
const double J_c5 = coordinate_dofs[1] * FE0_C0_D001_Q15[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D001_Q15[0][0][5 - 4] + coordinate_dofs[7] * FE0_C0_D001_Q15[0][0][6 - 4] + coordinate_dofs[10] * FE0_C0_D001_Q15[0][0][7 - 4];
const double J_c7 = coordinate_dofs[2] * FE0_C0_D010_Q15[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D010_Q15[0][0][9 - 8] + coordinate_dofs[8] * FE0_C0_D010_Q15[0][0][10 - 8];
const double J_c0 = coordinate_dofs[0] * FE0_C0_D100_Q15[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D100_Q15[0][0][1 - 0];
const double J_c1 = coordinate_dofs[0] * FE0_C0_D010_Q15[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D010_Q15[0][0][1 - 0] + coordinate_dofs[6] * FE0_C0_D010_Q15[0][0][2 - 0];
const double J_c6 = coordinate_dofs[2] * FE0_C0_D100_Q15[0][0][8 - 8] + coordinate_dofs[5] * FE0_C0_D100_Q15[0][0][9 - 8];
const double J_c3 = coordinate_dofs[1] * FE0_C0_D100_Q15[0][0][4 - 4] + coordinate_dofs[4] * FE0_C0_D100_Q15[0][0][5 - 4];
const double J_c2 = coordinate_dofs[0] * FE0_C0_D001_Q15[0][0][0 - 0] + coordinate_dofs[3] * FE0_C0_D001_Q15[0][0][1 - 0] + coordinate_dofs[6] * FE0_C0_D001_Q15[0][0][2 - 0] + coordinate_dofs[9] * FE0_C0_D001_Q15[0][0][3 - 0];
double sp15[84];
sp15[0] = J_c4 * J_c8;
sp15[1] = J_c5 * J_c7;
sp15[2] = -1 * sp15[1];
sp15[3] = sp15[0] + sp15[2];
sp15[4] = J_c0 * sp15[3];
sp15[5] = J_c5 * J_c6;
sp15[6] = J_c3 * J_c8;
sp15[7] = -1 * sp15[6];
sp15[8] = sp15[5] + sp15[7];
sp15[9] = J_c1 * sp15[8];
sp15[10] = sp15[4] + sp15[9];
sp15[11] = J_c3 * J_c7;
sp15[12] = J_c4 * J_c6;
sp15[13] = -1 * sp15[12];
sp15[14] = sp15[11] + sp15[13];
sp15[15] = J_c2 * sp15[14];
sp15[16] = sp15[10] + sp15[15];
sp15[17] = sp15[3] / sp15[16];
sp15[18] = -1 * J_c8;
sp15[19] = J_c3 * sp15[18];
sp15[20] = sp15[5] + sp15[19];
sp15[21] = sp15[20] / sp15[16];
sp15[22] = sp15[14] / sp15[16];
sp15[23] = sp15[22] * sp15[21];
sp15[24] = sp15[22] * sp15[22];
sp15[25] = sp15[22] * sp15[17];
sp15[26] = sp15[17] * sp15[21];
sp15[27] = sp15[17] * sp15[17];
sp15[28] = sp15[21] * sp15[21];
sp15[29] = J_c2 * J_c7;
sp15[30] = -1 * J_c1;
sp15[31] = J_c8 * sp15[30];
sp15[32] = sp15[29] + sp15[31];
sp15[33] = sp15[32] / sp15[16];
sp15[34] = J_c0 * J_c8;
sp15[35] = -1 * J_c2;
sp15[36] = J_c6 * sp15[35];
sp15[37] = sp15[34] + sp15[36];
sp15[38] = sp15[37] / sp15[16];
sp15[39] = J_c1 * J_c6;
sp15[40] = J_c0 * J_c7;
sp15[41] = -1 * sp15[40];
sp15[42] = sp15[39] + sp15[41];
sp15[43] = sp15[42] / sp15[16];
sp15[44] = sp15[38] * sp15[43];
sp15[45] = sp15[43] * sp15[43];
sp15[46] = sp15[33] * sp15[43];
sp15[47] = sp15[33] * sp15[38];
sp15[48] = sp15[33] * sp15[33];
sp15[49] = sp15[38] * sp15[38];
sp15[50] = sp15[48] + sp15[27];
sp15[51] = sp15[47] + sp15[26];
sp15[52] = sp15[46] + sp15[25];
sp15[53] = sp15[49] + sp15[28];
sp15[54] = sp15[44] + sp15[23];
sp15[55] = sp15[24] + sp15[45];
sp15[56] = J_c1 * J_c5;
sp15[57] = J_c2 * J_c4;
sp15[58] = -1 * sp15[57];
sp15[59] = sp15[56] + sp15[58];
sp15[60] = sp15[59] / sp15[16];
sp15[61] = J_c2 * J_c3;
sp15[62] = J_c0 * J_c5;
sp15[63] = -1 * sp15[62];
sp15[64] = sp15[61] + sp15[63];
sp15[65] = sp15[64] / sp15[16];
sp15[66] = J_c0 * J_c4;
sp15[67] = J_c1 * J_c3;
sp15[68] = -1 * sp15[67];
sp15[69] = sp15[66] + sp15[68];
sp15[70] = sp15[69] / sp15[16];
sp15[71] = sp15[70] * sp15[65];
sp15[72] = sp15[70] * sp15[70];
sp15[73] = sp15[70] * sp15[60];
sp15[74] = sp15[60] * sp15[65];
sp15[75] = sp15[60] * sp15[60];
sp15[76] = sp15[65] * sp15[65];
sp15[77] = sp15[50] + sp15[75];
sp15[78] = sp15[51] + sp15[74];
sp15[79] = sp15[52] + sp15[73];
sp15[80] = sp15[53] + sp15[76];
sp15[81] = sp15[54] + sp15[71];
sp15[82] = sp15[55] + sp15[72];
sp15[83] = fabs(sp15[16]);
for (int iq = 0; iq < 15; ++iq)
{
    // Quadrature loop body setup (num_points=15)
    // Section for geometrically varying computations
    double sv15[7];
    sv15[0] = weights15[iq] * sp15[83];
    sv15[1] = sp15[79] * sv15[0];
    sv15[2] = sp15[81] * sv15[0];
    sv15[3] = sp15[78] * sv15[0];
    sv15[4] = sp15[80] * sv15[0];
    sv15[5] = sp15[77] * sv15[0];
    sv15[6] = sp15[82] * sv15[0];
    for (int ia0 = 0; ia0 < 10; ++ia0)
    {
        for (int ia1 = 0; ia1 < 10; ++ia1)
        {
            A[ia0][ia1] += sv15[0] * FE1_C0_Q15[0][iq][ia0 - 0] * FE1_C0_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[5] * FE1_C0_D100_Q15[0][iq][ia0 - 0] * FE1_C0_D100_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[3] * FE1_C0_D100_Q15[0][iq][ia0 - 0] * FE1_C0_D010_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[1] * FE1_C0_D100_Q15[0][iq][ia0 - 0] * FE1_C0_D001_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[3] * FE1_C0_D010_Q15[0][iq][ia0 - 0] * FE1_C0_D100_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[4] * FE1_C0_D010_Q15[0][iq][ia0 - 0] * FE1_C0_D010_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[2] * FE1_C0_D010_Q15[0][iq][ia0 - 0] * FE1_C0_D001_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[1] * FE1_C0_D001_Q15[0][iq][ia0 - 0] * FE1_C0_D100_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[2] * FE1_C0_D001_Q15[0][iq][ia0 - 0] * FE1_C0_D010_Q15[0][iq][ia1 - 0];
            A[ia0][ia1] += sv15[6] * FE1_C0_D001_Q15[0][iq][ia0 - 0] * FE1_C0_D001_Q15[0][iq][ia1 - 0];
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
    double buffer_arg0_0[10][10]  = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec);
    MatSetValuesLocal(arg0_0_0, 10, arg0_0_map0_0 + i * 10,
                                             10, arg0_0_map1_0 + i * 10,
                                             (const PetscScalar *)buffer_arg0_0,
                                             ADD_VALUES);;
  }
}
        
        