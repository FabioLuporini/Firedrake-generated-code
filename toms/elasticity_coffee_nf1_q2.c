
        #include <petsc.h>
        #include <stdbool.h>
        #include <math.h>
        
        


        
            #include <immintrin.h>

            
            // This code is generated visiting a COFFEE AST

#include "firedrake_geometry.h"


static inline void form_cell_integral_0_otherwise (double A[30][30] , double** coordinate_dofs , double** w0 )
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
  
  
  static const double W14[14]  = {0.00317460317460317, 0.00317460317460317, 0.00317460317460317, 0.00317460317460317, 0.00317460317460317, 0.00317460317460317, 0.0147649707904968, 0.0147649707904968, 0.0147649707904968, 0.0147649707904968, 0.0221397911142651, 0.0221397911142651, 0.0221397911142651, 0.0221397911142651};
  static const double FE2_C1_D100[14][30]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2, 2, -2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 0.0, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2, 0.0, -2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 0.0, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 1.79367881729755, 0.0, 0.0, 0.0, 0.402107060900818, 0.402107060900818, -0.402107060900818, -0.402107060900818, -2.39157175639673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.79367881729755, -0.597892939099182, 0.0, 0.0, 0.0, 0.402107060900818, 0.402107060900818, -0.402107060900818, -0.402107060900818, 2.39157175639673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, -0.597892939099182, 0.0, 0.0, 0.0, 2.79367881729755, 0.402107060900818, -2.79367881729755, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099183, -0.597892939099185, 0.0, 0.0, 0.0, 0.402107060900818, 2.79367881729755, -0.402107060900818, -2.79367881729755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, -0.772474481918307, 0.0, 0.0, 0.0, 1.25749149397277, 1.25749149397277, -1.25749149397277, -1.25749149397277, 1.02996597589107, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.772474481918307, 0.257491493972768, 0.0, 0.0, 0.0, 1.25749149397277, 1.25749149397277, -1.25749149397277, -1.25749149397277, -1.02996597589108, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.257491493972768, 0.0, 0.0, 0.0, 0.227525518081694, 1.25749149397277, -0.227525518081694, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972769, 0.257491493972768, 0.0, 0.0, 0.0, 1.25749149397277, 0.227525518081694, -1.25749149397277, -0.227525518081694, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE2_C1_D010[14][30]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2, 0.0, 0.0, -2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 2, 0.0, 2, -2, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.00000000000001, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 2.00000000000001, 0.0, 2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 2, 0.0, 0.0, -2, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 0.0, -0.597892939099185, 0.0, 0.402107060900819, 0.0, 2.79367881729755, -0.402107060900818, 0.0, -2.79367881729755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.79367881729755, 0.0, -0.597892939099185, 0.0, 0.402107060900819, 0.0, 0.402107060900823, -0.402107060900818, 2.39157175639673, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099185, 0.0, -0.597892939099184, 0.0, 2.79367881729755, 0.0, 0.402107060900821, -2.79367881729755, 0.0, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 0.0, 1.79367881729754, 0.0, 0.402107060900819, 0.0, 0.402107060900825, -0.402107060900818, -2.39157175639673, -0.402107060900819, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.0, 0.257491493972766, 0.0, 1.25749149397277, 0.0, 0.227525518081699, -1.25749149397277, 0.0, -0.227525518081694, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.772474481918307, 0.0, 0.257491493972766, 0.0, 1.25749149397277, 0.0, 1.25749149397277, -1.25749149397277, -1.02996597589108, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972767, 0.0, 0.257491493972765, 0.0, 0.227525518081694, 0.0, 1.25749149397278, -0.227525518081694, 0.0, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972767, 0.0, -0.772474481918308, 0.0, 1.25749149397277, 0.0, 1.25749149397277, -1.25749149397277, 1.02996597589107, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE2_D010[14][30]  = {{1.0, 0.0, 1.0, 0.0, 2, 0.0, 0.0, -2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {1.0, 0.0, -1.0, 0.0, 2, 0.0, 2, -2, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.00000000000001, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 2.00000000000001, 0.0, 2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, 0.0, -1.0, 0.0, 2, 0.0, 0.0, -2, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099184, 0.0, -0.597892939099185, 0.0, 0.402107060900819, 0.0, 2.79367881729755, -0.402107060900818, 0.0, -2.79367881729755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.79367881729755, 0.0, -0.597892939099185, 0.0, 0.402107060900819, 0.0, 0.402107060900823, -0.402107060900818, 2.39157175639673, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099185, 0.0, -0.597892939099184, 0.0, 2.79367881729755, 0.0, 0.402107060900821, -2.79367881729755, 0.0, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099184, 0.0, 1.79367881729754, 0.0, 0.402107060900819, 0.0, 0.402107060900825, -0.402107060900818, -2.39157175639673, -0.402107060900819, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972768, 0.0, 0.257491493972766, 0.0, 1.25749149397277, 0.0, 0.227525518081699, -1.25749149397277, 0.0, -0.227525518081694, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.772474481918307, 0.0, 0.257491493972766, 0.0, 1.25749149397277, 0.0, 1.25749149397277, -1.25749149397277, -1.02996597589108, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972767, 0.0, 0.257491493972765, 0.0, 0.227525518081694, 0.0, 1.25749149397278, -0.227525518081694, 0.0, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972767, 0.0, -0.772474481918308, 0.0, 1.25749149397277, 0.0, 1.25749149397277, -1.25749149397277, 1.02996597589107, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE2_C2_D100[14][30]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2, 2, -2, -2, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 0.0, -2}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2, 0.0, -2, -2}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 2.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 0.0, 2}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 1.79367881729755, 0.0, 0.0, 0.0, 0.402107060900818, 0.402107060900818, -0.402107060900818, -0.402107060900818, -2.39157175639673}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.79367881729755, -0.597892939099182, 0.0, 0.0, 0.0, 0.402107060900818, 0.402107060900818, -0.402107060900818, -0.402107060900818, 2.39157175639673}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, -0.597892939099182, 0.0, 0.0, 0.0, 2.79367881729755, 0.402107060900818, -2.79367881729755, -0.402107060900818, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099183, -0.597892939099185, 0.0, 0.0, 0.0, 0.402107060900818, 2.79367881729755, -0.402107060900818, -2.79367881729755, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, -0.772474481918307, 0.0, 0.0, 0.0, 1.25749149397277, 1.25749149397277, -1.25749149397277, -1.25749149397277, 1.02996597589107}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.772474481918307, 0.257491493972768, 0.0, 0.0, 0.0, 1.25749149397277, 1.25749149397277, -1.25749149397277, -1.25749149397277, -1.02996597589108}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.257491493972768, 0.0, 0.0, 0.0, 0.227525518081694, 1.25749149397277, -0.227525518081694, -1.25749149397277, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972769, 0.257491493972768, 0.0, 0.0, 0.0, 1.25749149397277, 0.227525518081694, -1.25749149397277, -0.227525518081694, 0.0}};
  static const double FE2_D100[14][30]  = {{1.0, -1.0, 0.0, 0.0, 0.0, 2, 2, -2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {1.0, 1.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 0.0, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2, 0.0, -2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, -1.0, 0.0, 0.0, 0.0, 2, 0.0, -2, 0.0, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099184, 1.79367881729755, 0.0, 0.0, 0.0, 0.402107060900818, 0.402107060900818, -0.402107060900818, -0.402107060900818, -2.39157175639673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.79367881729755, -0.597892939099182, 0.0, 0.0, 0.0, 0.402107060900818, 0.402107060900818, -0.402107060900818, -0.402107060900818, 2.39157175639673, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099184, -0.597892939099182, 0.0, 0.0, 0.0, 2.79367881729755, 0.402107060900818, -2.79367881729755, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099183, -0.597892939099185, 0.0, 0.0, 0.0, 0.402107060900818, 2.79367881729755, -0.402107060900818, -2.79367881729755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972768, -0.772474481918307, 0.0, 0.0, 0.0, 1.25749149397277, 1.25749149397277, -1.25749149397277, -1.25749149397277, 1.02996597589107, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.772474481918307, 0.257491493972768, 0.0, 0.0, 0.0, 1.25749149397277, 1.25749149397277, -1.25749149397277, -1.25749149397277, -1.02996597589108, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972768, 0.257491493972768, 0.0, 0.0, 0.0, 0.227525518081694, 1.25749149397277, -0.227525518081694, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972769, 0.257491493972768, 0.0, 0.0, 0.0, 1.25749149397277, 0.227525518081694, -1.25749149397277, -0.227525518081694, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
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
  static const double FE2_C2_D010[14][30]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2, 0.0, 0.0, -2, -2, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 2, 0.0, 2, -2, 0.0, -2.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.00000000000001, 0.0, -2.0, -2.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 2.00000000000001, 0.0, 2, -2}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 2, 0.0, 0.0, -2, 2.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 0.0, -0.597892939099185, 0.0, 0.402107060900819, 0.0, 2.79367881729755, -0.402107060900818, 0.0, -2.79367881729755}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.79367881729755, 0.0, -0.597892939099185, 0.0, 0.402107060900819, 0.0, 0.402107060900823, -0.402107060900818, 2.39157175639673, -0.402107060900818}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099185, 0.0, -0.597892939099184, 0.0, 2.79367881729755, 0.0, 0.402107060900821, -2.79367881729755, 0.0, -0.402107060900818}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 0.0, 1.79367881729754, 0.0, 0.402107060900819, 0.0, 0.402107060900825, -0.402107060900818, -2.39157175639673, -0.402107060900819}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.0, 0.257491493972766, 0.0, 1.25749149397277, 0.0, 0.227525518081699, -1.25749149397277, 0.0, -0.227525518081694}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.772474481918307, 0.0, 0.257491493972766, 0.0, 1.25749149397277, 0.0, 1.25749149397277, -1.25749149397277, -1.02996597589108, -1.25749149397277}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972767, 0.0, 0.257491493972765, 0.0, 0.227525518081694, 0.0, 1.25749149397278, -0.227525518081694, 0.0, -1.25749149397277}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972767, 0.0, -0.772474481918308, 0.0, 1.25749149397277, 0.0, 1.25749149397277, -1.25749149397277, 1.02996597589107, -1.25749149397277}};
  static const double FE2_C2_D001[14][30]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, -2, -2.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.00000000000001, 0.0, -2, 0.0, -2.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 2, 2.00000000000001, 0.0, 0.0, -2.0, -2}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 2.00000000000001, 0.0, 2.0, 0.0, -2}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 2, 0.0, 0.0, 2, -2, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 0.0, 0.0, -0.597892939099186, 0.40210706090082, 2.79367881729756, 0.0, 0.0, -0.402107060900818, -2.79367881729755}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.79367881729755, 0.0, 0.0, -0.597892939099187, 0.40210706090082, 0.402107060900823, 0.0, 2.39157175639673, -0.402107060900818, -0.402107060900818}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 0.0, 0.0, 1.79367881729754, 0.402107060900815, 0.402107060900826, 0.0, -2.39157175639673, -0.402107060900816, -0.402107060900819}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099183, 0.0, 0.0, -0.597892939099186, 2.79367881729755, 0.402107060900822, 0.0, 0.0, -2.79367881729755, -0.402107060900818}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.0, 0.0, 0.257491493972765, 1.25749149397277, 0.227525518081701, 0.0, 0.0, -1.25749149397277, -0.227525518081694}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.772474481918307, 0.0, 0.0, 0.257491493972765, 1.25749149397277, 1.25749149397278, 0.0, -1.02996597589108, -1.25749149397277, -1.25749149397277}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.0, 0.0, -0.772474481918311, 1.25749149397277, 1.25749149397278, 0.0, 1.02996597589107, -1.25749149397277, -1.25749149397277}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.0, 0.0, 0.257491493972765, 0.227525518081694, 1.25749149397278, 0.0, 0.0, -0.227525518081694, -1.25749149397277}};
  static const double FE2_D001[14][30]  = {{1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, -2, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {1.0, 0.0, 0.0, 1.0, 0.0, 2.00000000000001, 0.0, -2, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {1.0, 0.0, 0.0, -1.0, 2, 2.00000000000001, 0.0, 0.0, -2.0, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, 0.0, 0.0, -1.0, 0.0, 2.00000000000001, 0.0, 2.0, 0.0, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, 0.0, 0.0, -1.0, 2, 0.0, 0.0, 2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099184, 0.0, 0.0, -0.597892939099186, 0.40210706090082, 2.79367881729756, 0.0, 0.0, -0.402107060900818, -2.79367881729755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-1.79367881729755, 0.0, 0.0, -0.597892939099187, 0.40210706090082, 0.402107060900823, 0.0, 2.39157175639673, -0.402107060900818, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099184, 0.0, 0.0, 1.79367881729754, 0.402107060900815, 0.402107060900826, 0.0, -2.39157175639673, -0.402107060900816, -0.402107060900819, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.597892939099183, 0.0, 0.0, -0.597892939099186, 2.79367881729755, 0.402107060900822, 0.0, 0.0, -2.79367881729755, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972768, 0.0, 0.0, 0.257491493972765, 1.25749149397277, 0.227525518081701, 0.0, 0.0, -1.25749149397277, -0.227525518081694, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.772474481918307, 0.0, 0.0, 0.257491493972765, 1.25749149397277, 1.25749149397278, 0.0, -1.02996597589108, -1.25749149397277, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972768, 0.0, 0.0, -0.772474481918311, 1.25749149397277, 1.25749149397278, 0.0, 1.02996597589107, -1.25749149397277, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {-0.257491493972768, 0.0, 0.0, 0.257491493972765, 0.227525518081694, 1.25749149397278, 0.0, 0.0, -0.227525518081694, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  static const double FE2_C1_D001[14][30]  = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, -2, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 2.00000000000001, 0.0, -2, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 2, 2.00000000000001, 0.0, 0.0, -2.0, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 2.00000000000001, 0.0, 2.0, 0.0, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 2, 0.0, 0.0, 2, -2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 0.0, 0.0, -0.597892939099186, 0.40210706090082, 2.79367881729756, 0.0, 0.0, -0.402107060900818, -2.79367881729755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.79367881729755, 0.0, 0.0, -0.597892939099187, 0.40210706090082, 0.402107060900823, 0.0, 2.39157175639673, -0.402107060900818, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099184, 0.0, 0.0, 1.79367881729754, 0.402107060900815, 0.402107060900826, 0.0, -2.39157175639673, -0.402107060900816, -0.402107060900819, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.597892939099183, 0.0, 0.0, -0.597892939099186, 2.79367881729755, 0.402107060900822, 0.0, 0.0, -2.79367881729755, -0.402107060900818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.0, 0.0, 0.257491493972765, 1.25749149397277, 0.227525518081701, 0.0, 0.0, -1.25749149397277, -0.227525518081694, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.772474481918307, 0.0, 0.0, 0.257491493972765, 1.25749149397277, 1.25749149397278, 0.0, -1.02996597589108, -1.25749149397277, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.0, 0.0, -0.772474481918311, 1.25749149397277, 1.25749149397278, 0.0, 1.02996597589107, -1.25749149397277, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.257491493972768, 0.0, 0.0, 0.257491493972765, 0.227525518081694, 1.25749149397278, 0.0, 0.0, -0.227525518081694, -1.25749149397277, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  
  
  for (int ip  = 0; ip < 14; ip += 1)
  {
    double F0  = 0.0;
    
    for (int r  = 0; r < 10; r += 1)
    {
      F0 += (w0[r][0] * FE0[ip][r]);
      
    }
    double IP_1_1_0  = (det * W14[ip] * F0 * 0.25);
    double IP_J_1_3_0[30]  = {0.0};
    double IP_J_1_3_1[30]  = {0.0};
    double IP_J_1_3_2[30]  = {0.0};
    double IP_J_1_3_3[30]  = {0.0};
    double IP_J_1_3_4[30]  = {0.0};
    double IP_J_1_3_5[30]  = {0.0};
    double IP_K_1_3_0[30]  = {0.0};
    double IP_K_1_3_1[30]  = {0.0};
    double IP_K_1_3_2[30]  = {0.0};
    double IP_K_1_3_3[30]  = {0.0};
    double IP_K_1_3_4[30]  = {0.0};
    double IP_K_1_3_5[30]  = {0.0};
    
    for (int k  = 0; k < 10; k += 1)
    {
      IP_J_1_3_0[k+10] = (K[1] * FE2_C1_D100[ip][k+10]) + (K[4] * FE2_C1_D010[ip][k+10]) + (K[7] * FE2_C1_D001[ip][k+10]);
      IP_J_1_3_1[k] = (K[2] * FE2_D100[ip][k]) + (K[5] * FE2_D010[ip][k]) + (K[8] * FE2_D001[ip][k]) + (K[0] * FE2_C2_D100[ip][k]) + (K[3] * FE2_C2_D010[ip][k]) + (K[6] * FE2_C2_D001[ip][k]);
      IP_J_1_3_1[k+20] = (K[2] * FE2_D100[ip][k+20]) + (K[5] * FE2_D010[ip][k+20]) + (K[8] * FE2_D001[ip][k+20]) + (K[0] * FE2_C2_D100[ip][k+20]) + (K[3] * FE2_C2_D010[ip][k+20]) + (K[6] * FE2_C2_D001[ip][k+20]);
      IP_J_1_3_4[k] = (K[0] * FE2_D100[ip][k]) + (K[3] * FE2_D010[ip][k]) + (K[6] * FE2_D001[ip][k]);
      IP_J_1_3_5[k+20] = (K[2] * FE2_C2_D100[ip][k+20]) + (K[5] * FE2_C2_D010[ip][k+20]) + (K[8] * FE2_C2_D001[ip][k+20]);
      IP_K_1_3_0[k] = ((FE2_D100[ip][k] * 2 * K[0]) + (FE2_D010[ip][k] * 2 * K[3]) + (FE2_D001[ip][k] * 2 * K[6])) * 2 * IP_1_1_0;
      IP_K_1_3_1[k+20] = ((FE2_C2_D100[ip][k+20] * 2 * K[2]) + (FE2_C2_D010[ip][k+20] * 2 * K[5]) + (FE2_C2_D001[ip][k+20] * 2 * K[8])) * 2 * IP_1_1_0;
      IP_K_1_3_2[k+20] = ((FE2_D100[ip][k+20] * K[2]) + (FE2_D010[ip][k+20] * K[5]) + (FE2_D001[ip][k+20] * K[8]) + (FE2_C2_D100[ip][k+20] * K[0]) + (FE2_C2_D010[ip][k+20] * K[3]) + (FE2_C2_D001[ip][k+20] * K[6])) * 2 * IP_1_1_0;
      IP_K_1_3_2[k] = ((FE2_D100[ip][k] * K[2]) + (FE2_D010[ip][k] * K[5]) + (FE2_D001[ip][k] * K[8]) + (FE2_C2_D100[ip][k] * K[0]) + (FE2_C2_D010[ip][k] * K[3]) + (FE2_C2_D001[ip][k] * K[6])) * 2 * IP_1_1_0;
      IP_K_1_3_3[k+10] = ((FE2_C1_D100[ip][k+10] * 2 * K[1]) + (FE2_C1_D010[ip][k+10] * 2 * K[4]) + (FE2_C1_D001[ip][k+10] * 2 * K[7])) * 2 * IP_1_1_0;
      
    }
    
    for (int k  = 0; k < 20; k += 1)
    {
      IP_J_1_3_2[k] = (K[1] * FE2_D100[ip][k]) + (K[4] * FE2_D010[ip][k]) + (K[7] * FE2_D001[ip][k]) + (K[0] * FE2_C1_D100[ip][k]) + (K[3] * FE2_C1_D010[ip][k]) + (K[6] * FE2_C1_D001[ip][k]);
      IP_J_1_3_3[k+10] = (K[1] * FE2_C2_D100[ip][k+10]) + (K[4] * FE2_C2_D010[ip][k+10]) + (K[7] * FE2_C2_D001[ip][k+10]) + (K[2] * FE2_C1_D100[ip][k+10]) + (K[5] * FE2_C1_D010[ip][k+10]) + (K[8] * FE2_C1_D001[ip][k+10]);
      IP_K_1_3_4[k+10] = ((FE2_C2_D100[ip][k+10] * K[1]) + (FE2_C2_D010[ip][k+10] * K[4]) + (FE2_C2_D001[ip][k+10] * K[7]) + (FE2_C1_D100[ip][k+10] * K[2]) + (FE2_C1_D010[ip][k+10] * K[5]) + (FE2_C1_D001[ip][k+10] * K[8])) * 2 * IP_1_1_0;
      IP_K_1_3_5[k] = ((FE2_D100[ip][k] * K[1]) + (FE2_D010[ip][k] * K[4]) + (FE2_D001[ip][k] * K[7]) + (FE2_C1_D100[ip][k] * K[0]) + (FE2_C1_D010[ip][k] * K[3]) + (FE2_C1_D001[ip][k] * K[6])) * 2 * IP_1_1_0;
      
    }
    
    for (int j  = 0; j < 10; j += 1)
    {
      
      for (int k  = 0; k < 10; k += 1)
      {
        A[j+10][k+10] += ((IP_K_1_3_3[k+10] * IP_J_1_3_0[j+10]));
        A[j+20][k+20] += ((IP_K_1_3_1[k+20] * IP_J_1_3_5[j+20])) + ((IP_K_1_3_2[k+20] * IP_J_1_3_1[j+20]));
        A[j][k] += ((IP_K_1_3_2[k] * IP_J_1_3_1[j])) + ((IP_K_1_3_0[k] * IP_J_1_3_4[j]));
        A[j+20][k] += ((IP_K_1_3_2[k] * IP_J_1_3_1[j+20]));
        A[j][k+20] += ((IP_K_1_3_2[k+20] * IP_J_1_3_1[j]));
        
      }
      
    }
    
    for (int j  = 0; j < 20; j += 1)
    {
      
      for (int k  = 0; k < 20; k += 1)
      {
        A[j+10][k+10] += ((IP_K_1_3_4[k+10] * IP_J_1_3_3[j+10]));
        A[j][k] += ((IP_K_1_3_5[k] * IP_J_1_3_2[j]));
        
      }
      
    }
    
  }
  
}
            

        
        void wrap_form_cell_integral_0_otherwise(int start, int end,
                      Mat arg0_0_, int *arg0_0_map0_0, int *arg0_0_map1_0, double *arg1_0, int *arg1_0_map0_0, double *arg2_0, int *arg2_0_map0_0
                      ) {
  Mat arg0_0_0 = arg0_0_;
  double *arg1_0_vec[12];
    double *arg2_0_vec[10];
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
    arg2_0_vec[0] = arg2_0 + (arg2_0_map0_0[i * 10 + 0])* 1;
    arg2_0_vec[1] = arg2_0 + (arg2_0_map0_0[i * 10 + 1])* 1;
    arg2_0_vec[2] = arg2_0 + (arg2_0_map0_0[i * 10 + 2])* 1;
    arg2_0_vec[3] = arg2_0 + (arg2_0_map0_0[i * 10 + 3])* 1;
    arg2_0_vec[4] = arg2_0 + (arg2_0_map0_0[i * 10 + 4])* 1;
    arg2_0_vec[5] = arg2_0 + (arg2_0_map0_0[i * 10 + 5])* 1;
    arg2_0_vec[6] = arg2_0 + (arg2_0_map0_0[i * 10 + 6])* 1;
    arg2_0_vec[7] = arg2_0 + (arg2_0_map0_0[i * 10 + 7])* 1;
    arg2_0_vec[8] = arg2_0 + (arg2_0_map0_0[i * 10 + 8])* 1;
    arg2_0_vec[9] = arg2_0 + (arg2_0_map0_0[i * 10 + 9])* 1;
    double buffer_arg0_0[30][30]  = {{0.0}};
    form_cell_integral_0_otherwise(buffer_arg0_0, arg1_0_vec, arg2_0_vec);
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
        
        