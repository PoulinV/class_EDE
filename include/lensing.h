/** @file lensing.h Documented includes for spectra module */

#ifndef __LENSING__
#define __LENSING__

#include "spectra.h"

/**
 * Structure containing everything about lensed spectra that other modules need to know.
 *
 * Once initialized by lensing_init(), contains a table of all lensed
 * C_l's for the all modes (scalar/tensor), all types (TT, TE...),
 * and all pairs of initial conditions (adiabatic, isocurvatures...).
 * FOR THE MOMENT, ASSUME ONLY SCALAR & ADIABATIC
 */

struct lensing {

  /** @name - input parameters initialized by user in input module
   *  (all other quantitites are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */
  
  //@{

  short has_lensed_cls; /**< do we need to compute lensed Cl's at all ? */

  //@}

  /** @name - information on number of type of C_l's (TT, TE...) */

  //@{

  int has_tt; /**< do we want lensed C_l^TT ? (T = temperature) */
  int has_ee; /**< do we want lensed C_l^EE ? (E = E-polarization) */
  int has_te; /**< do we want lensed C_l^TE ? */
  int has_bb; /**< do we want C_l^BB ? (B = B-polarization) */
  int has_pp; /**< do we want C_l^phi-phi ? (phi = CMB lensing potential) */
  int has_tp; /**< do we want C_l^T-phi ? */

  int index_lt_tt; /**< index for type C_l^TT */
  int index_lt_ee; /**< index for type C_l^EE */
  int index_lt_te; /**< index for type C_l^TE */
  int index_lt_bb; /**< index for type C_l^BB */
  int index_lt_pp; /**< index for type C_l^phi-phi */
  int index_lt_tp; /**< index for type C_l^T-phi */

  int lt_size; /**< number of C_l types requested */

  //@}

  /** @name - table of pre-computed C_l values, and related quantitites */

  //@{

  int l_unlensed_max;    /**< last multipole in all calculations (same as in spectra module)*/

  int l_lensed_max;    /**< last multipole at which lensed spactra are computed */

  double * cl_lensed;   /**< table of anisotropy spectra for each
			   multipole and types, 
			   cl[index_l * ple->lt_size + index_lt]; 
			   index_l=l goes from 0 to l_lensed_max, but for index_l=0,1
			   cl_lensed is set to zero; its size is
			   (l_max+1)*psp->ct_size */


  //@}

  /** @name - technical parameters */

  //@{

  short lensing_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif

  int lensing_cl_at_l(
                      struct lensing * ple,
                      int l,
                      double * cl_lensed
                      );

  int lensing_init(
		   struct precision * ppr,
                   struct perturbs * ppt,
                   struct spectra * psp,
                   struct lensing * ple
                   );

  int lensing_free(
                   struct lensing * ple
                   );
    
  int lensing_indices(
		      struct precision * ppr,
                      struct perturbs * ppt,
                      struct spectra * psp,
                      struct lensing * ple
                      );
  
  int lensing_lensed_cl_tt(
                        double *ksi, 
                        double **d00,
                        double *w8,
                        int nmu,
                        struct lensing * ple
                        );
  
  int lensing_lensed_cl_te(
                           double *ksiX, 
                           double **d20,
                           double *w8,
                           int nmu,
                           struct lensing * ple
                           );
  
  int lensing_X000(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X000
                   );
    
  int lensing_Xp000(
                    double * mu,
                    int num_mu,
                    int lmax,
                    double * sigma2,
                    double ** Xp000
                    );
  
  int lensing_X220(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X220
                   );
  
  int lensing_X022(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X022
                   );

  int lensing_Xp022(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** Xp022
                   );
  
  int lensing_X121(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X121
                   );
  
  int lensing_X132(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X132
                   );
  
  int lensing_X242(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double * sigma2,
                   double ** X242
                   );
  
  int lensing_d00(
                  double * mu,
                  int num_mu,
                  int lmax,
                  double ** d00
                  );

  int lensing_d11(
                  double * mu,
                  int num_mu,
                  int lmax,
                  double ** d11
                  );

  int lensing_d1m1(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d1m1
                   );

  int lensing_d2m2(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d2m2
                   );
  
  int lensing_d22(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d22
                   );

  int lensing_d20(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d20
                   );
  
  int lensing_d31(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d3m1
                   );
  
  int lensing_d3m1(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d3m1
                   );
  
  int lensing_d3m3(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d3m3
                   );
  
  int lensing_d40(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d40
                   );
  
  int lensing_d4m2(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d4m2
                   );
  
  int lensing_d4m4(
                   double * mu,
                   int num_mu,
                   int lmax,
                   double ** d4m4
                   );
  
  int lensing_gauss_legendre(
                             double *mu,
                             double *w8,
                             int nmu
                             );
    
    
#ifdef __cplusplus
}
#endif

#endif
