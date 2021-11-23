/** @file lensing.h Documented includes for harmonic module */

#ifndef __LENSING__
#define __LENSING__

#include "harmonic.h"

/**
 * Structure containing everything about lensed spectra that other modules need to know.
 *
 * Once initialized by lensing_init(), contains a table of all lensed
 * \f$ C_l\f$'s for the all modes (scalar/tensor), all types (TT, TE...),
 * and all pairs of initial conditions (adiabatic, isocurvatures...).
 * FOR THE MOMENT, ASSUME ONLY SCALAR & ADIABATIC
 */
 
 /** DLM: enum defining whether cmb spectra should be calculated or else given externally */
enum cmb_spectra_type {
	internal_cmb,
	external_cmb
};
/** DLM: enum defining whether lensed cmb spectra should be calculated or else given externally */
enum cmb_spectra_lensed_type {
	internal_cmb_lensed,
	external_cmb_lensed
};
/** DLM: enum defining how the lensing temparature noise spectra should be computed */
enum temperature_noise_type {
	idealized_tn,
	external_tn
};
/** DLM: enum defining how the lensing polariation noise spectra should be computed */
enum polarization_noise_type {
	idealized_pn,
	external_pn
};
/** DLM: enum defining how the lensing polariation noise spectra should be computed */
enum lens_rec_noise_type {
	internal_rn,
	external_rn
};
/** DLM: enum defining whether the noise calculated carry units or not (default:unifull)*/
enum noise_type {
	unitfull,
	unitless
};
enum convergence_type{
	every,
	total
};
enum derv_type{
	lensed,
	delensed
};
/*--------------------------------*/

struct lensing {

  /** @name - input parameters initialized by user in input module
   *  (all other quantities are computed in this module, given these
   *  parameters and the content of the 'precision', 'background' and
   *  'thermodynamics' structures) */

  //@{
	/*----------------DLM-------------*/
	/*--------------------------------*/
	enum temperature_noise_type temperature_noise_type; /* DLM */
	enum polarization_noise_type polarization_noise_type; /* DLM */
	enum lens_rec_noise_type lens_rec_noise_type; /* DLM */
	enum cmb_spectra_type cmb_spectra_type; /* DLM */
	enum cmb_spectra_lensed_type cmb_spectra_lensed_type; /* DLM */

	enum convergence_type convergence_type; /* DLM */
	enum derv_type derv_type; /* DLM */

	enum noise_type noise_type; /* DLM */
	/*--------------------------------*/

	short has_lensed_cls; /**< do we need to compute lensed \f$ C_l\f$'s at all ? */

	short calculate_pderivaties;
	//@}

  /** @name - information on number of type of C_l's (TT, TE...) */

  //@{

	char* command_for_temp_noise_spec; /* DLM */
	char* command_for_polarization_noise_spec; /* DLM */
	char* command_for_lens_recon_noise_spec; /* DLM */

    char* command_for_external_cmb_spectra; /* DLM */
    char* command_for_external_lensed_cmb_spectra; /* DLM */

	int has_tt; /**< do we want lensed \f$ C_l^{TT}\f$? (T = temperature) */
	int has_ee; /**< do we want lensed \f$ C_l^{EE}\f$? (E = E-polarization) */
	int has_te; /**< do we want lensed \f$ C_l^{TE}\f$? */
	int has_bb; /**< do we want \f$ C_l^{BB}\f$? (B = B-polarization) */
	int has_pp; /**< do we want \f$ C_l^{\phi\phi}\f$? (\f$ \phi \f$ = CMB lensing potential) */
	int has_tp; /**< do we want \f$ C_l^{T\phi}\f$? */
	int has_dd; /**< do we want \f$ C_l^{dd}\f$? (d = matter density) */
	int has_td; /**< do we want \f$ C_l^{Td}\f$? */
	int has_ll; /**< do we want \f$ C_l^{ll}\f$? (l = lensing potential) */
	int has_tl; /**< do we want \f$ C_l^{Tl}\f$? */


	int index_lt_tt; /**< index for type \f$ C_l^{TT} \f$*/
	int index_lt_ee; /**< index for type \f$ C_l^{EE} \f$*/
	int index_lt_te; /**< index for type \f$ C_l^{TE} \f$*/
	int index_lt_bb; /**< index for type \f$ C_l^{BB} \f$*/
	int index_lt_pp; /**< index for type \f$ C_l^{\phi\phi} \f$*/
	int index_lt_tp; /**< index for type \f$ C_l^{T\phi} \f$*/
	int index_lt_dd; /**< index for type \f$ C_l^{dd} \f$*/
	int index_lt_td; /**< index for type \f$ C_l^{Td} \f$*/
	int index_lt_ll; /**< index for type \f$ C_l^{dd} \f$*/
	int index_lt_tl; /**< index for type \f$ C_l^{Td} \f$*/


	int lt_size; /**< number of \f$ C_l\f$ types requested */

	//@}

  /** @name - table of pre-computed C_l values, and related quantities */

  //@{

  int l_unlensed_max;    /**< last multipole in all calculations (same as in harmonic module)*/

  int l_lensed_max;    /**< last multipole at which lensed spectra are computed */

  /* interpolable version: */

  int l_size;       /**< number of l values */

  int * l_max_lt;    /**< last multipole (given as an input) at which
		    we want to output \f$ C_l \f$'s for a given mode and type */

  double * l;       /**< table of multipole values l[index_l] */
  double * l_unlensed;
  double * cl_lens; /**< table of anisotropy spectra for each
			   multipole and types,
			   cl[index_l * ple->lt_size + index_lt] */

  double ** cl_lens_derv_TT_TT;
  double ** cl_lens_derv_TE_TE;
  double ** cl_lens_derv_EE_EE;
  double ** cl_lens_derv_EE_BB;
  double ** cl_lens_derv_BB_EE;
  double ** cl_lens_derv_BB_BB;
    
  double ** ddcl_lens_derv_TT_TT;
  double ** ddcl_lens_derv_TE_TE;
  double ** ddcl_lens_derv_EE_EE;
  double ** ddcl_lens_derv_EE_BB;
  double ** ddcl_lens_derv_BB_EE;
  double ** ddcl_lens_derv_BB_BB;

  double ** cl_delens_derv_TT_TT;
  double ** cl_delens_derv_TE_TE;
  double ** cl_delens_derv_EE_EE;
  double ** cl_delens_derv_EE_BB;
  double ** cl_delens_derv_BB_EE;
  double ** cl_delens_derv_BB_BB;
    
  double ** ddcl_delens_derv_TT_TT;
  double ** ddcl_delens_derv_TE_TE;
  double ** ddcl_delens_derv_EE_EE;
  double ** ddcl_delens_derv_EE_BB;
  double ** ddcl_delens_derv_BB_EE;
  double ** ddcl_delens_derv_BB_BB;

  double * ddcl_lens; /**< second derivatives for interpolation */

  int dl_size;       /**< DLM: number of l values for the delensed spectra */

	short has_delensed_cls; /**< DLM: flag is the delensed spectra is to be calculated */

	short has_itr_delensing; /**< DLM: flag is the delensed spectra is to be calculated iteratively */

	int max_itr_steps;  /**< DLM: integer equal to the maximum steps for lensing noise reconstruction iteration */
	double convergence_criterion_itr;/**< DLM: double equal to the minimum desired percent ratio of the difference between iterations */

	short has_lens_noise_rcn; /**< DLM: flag if the lensing noise is available (otherwise will be calculated) */

	double * l_dl;       /**< DLM: table of multipole values l[index_l] for the delensed spectra */
	double * cl_delens;   /**< DLM: delensed spectra */
	double * ddcl_delens; /**< DLM: second derivatives for interpolation */

	double ** cl_dl_tt_pderv; /*DLM: and dl_size*dl_size matrix of cl spectra derivatives w.r.t. lensing spectra.*/
	double ** cl_dl_te_pderv; /*DLM: and dl_size*dl_size matrix of cl spectra derivatives w.r.t. lensing spectra.*/
	double ** cl_dl_ee_pderv; /*DLM: and dl_size*dl_size matrix of cl spectra derivatives w.r.t. lensing spectra.*/
	double ** cl_dl_bb_pderv; /*DLM: and dl_size*dl_size matrix of cl spectra derivatives w.r.t. lensing spectra.*/

	double ** ddcl_dl_tt_pderv; /*DLM: and dl_size*dl_size matrix of cl spectra derivatives w.r.t. lensing spectra.*/
	double ** ddcl_dl_te_pderv; /*DLM: and dl_size*dl_size matrix of cl spectra derivatives w.r.t. lensing spectra.*/
	double ** ddcl_dl_ee_pderv; /*DLM: and dl_size*dl_size matrix of cl spectra derivatives w.r.t. lensing spectra.*/
	double ** ddcl_dl_bb_pderv; /*DLM: and dl_size*dl_size matrix of cl spectra derivatives w.r.t. lensing spectra.*/

	int l_delensed_max; /* DLM */

	int dlt_size;  /**< DLM: number of delensed \f$ C_l\f$ types requested */

	int index_lt_dl_tt; /**< DLM: index for type \f$ C_l^{TT,delensed} \f$*/
	int index_lt_dl_ee; /**< DLM: index for type \f$ C_l^{EE,delensed} \f$*/
	int index_lt_dl_te; /**< DLM: index for type \f$ C_l^{TE,delensed} \f$*/
	int index_lt_dl_bb; /**< DLM: index for type \f$ C_l^{BB,delensed} \f$*/
	int index_lt_dl_pp; /**< DLM: index for type \f$ C_l^{\phi\phi} \f$*/
	int index_lt_dl_tp; /**< DLM: index for type \f$ C_l^{T\phi,delensed} \f$*/
	int index_lt_dl_dd; /**< DLM: index for type \f$ C_l^{dd} \f$*/
	int index_lt_dl_td; /**< DLM: index for type \f$ C_l^{Td,delensed} \f$*/
	int index_lt_dl_ll; /**< DLM: index for type \f$ C_l^{dd,delensed} \f$*/
	int index_lt_dl_tl; /**< DLM: index for type \f$ C_l^{Td,delensed} \f$*/

	int has_dl_tt; /**< DLM: TRUE if \f$ C_l^{TT,delensed} \f$*/
	int has_dl_ee; /**< DLM: TRUE if \f$ C_l^{EE,delensed} \f$*/
	int has_dl_te; /**< DLM: TRUE if \f$ C_l^{TE,delensed} \f$*/
	int has_dl_bb; /**< DLM: TRUE if \f$ C_l^{BB,delensed} \f$*/
	int has_dl_pp; /**< DLM: TRUE if \f$ C_l^{\phi\phi,delensed} \f$*/
	int has_dl_tp; /**< DLM: TRUE if \f$ C_l^{T\phi,delensed} \f$*/
	int has_dl_dd; /**< DLM: TRUE if \f$ C_l^{dd,delensed} \f$*/
	int has_dl_td; /**< DLM: TRUE if \f$ C_l^{Td,delensed} \f$*/
	int has_dl_ll; /**< DLM: TRUE if \f$ C_l^{dd,delensed} \f$*/
	int has_dl_tl; /**< DLM: TRUE if \f$ C_l^{Td,delensed} \f$*/

	short output_spectra_noise; /**< DLM: TRUE if the spectra noise is to be written */

	short output_derivatives; /**< DLM: TRUE if the derivatives of the spectra w.r.t. lensing spectrum is to be written to files*/

	int derv_binedges; /* DLM */
	
    int recon_mask_lmin_T; /* DLM: Used for specifying l_mask for lensing reconstruction */
    int recon_mask_lmax_T;
    int recon_mask_lmin_E;
    int recon_mask_lmax_E;
    int recon_mask_lmin_B;
    int recon_mask_lmax_B;

	double * l_rcn_ext, * pk_rcn_ext; /**< DLM: external lensing noise \f$*/

	double * l_tn, * pk_tn; /**< DLM: Temperature noise \f$*/
	double * l_pn, * pk_pn; /**< DLM: Polarization noise \f$*/

	short has_nl_all; /* DLM: TRUE if all spectra are used to generate the min-var*/
	short has_nl_diag; /* DLM: TRUE if only the varriances are used to generate the min-var*/
	short has_nl_eb; /* DLM: TRUE if only the varriances are used to generate the min-var*/

	short has_nl_all_itr; /* DLM: TRUE if all spectra are used to generate the min-var*/
	short has_nl_diag_itr; /* DLM: TRUE if only the varriances are used to generate the min-var*/
	short has_nl_eb_itr; /* DLM: TRUE if only the varriances are used to generate the min-var*/
	short has_nl_altr_itr; /* DLM: TRUE if only the varriances are used to generate the min-var*/

	double * nl_rcn; /* DLM */
	double * ddnl_rcn; /* DLM */

	short has_ln_rcn_tt; /**< DLM: True if have \f$ N_l^{\phi\phi,TT} \f$*/
	short has_ln_rcn_ee; /**< DLM: True if have \f$ N_l^{\phi\phi,EE} \f$*/
	short has_ln_rcn_te; /**< DLM: True if have \f$ N_l^{\phi\phi,TE} \f$*/
	short has_ln_rcn_bb; /**< DLM: True if have \f$ N_l^{\phi\phi,BB} \f$*/
	short has_ln_rcn_tb; /**< DLM: True if have \f$ N_l^{\phi\phi,TB} \f$*/
	short has_ln_rcn_eb; /**< DLM: True if have \f$ N_l^{\phi\phi,EB} \f$*/
	short has_ln_rcn_ttte;/**< DLM: True if have \f$ C_l^{\phi\phi,TTTE} \f$*/
	short has_ln_rcn_ttee;/**< DLM: True if have \f$ C_l^{\phi\phi,TTTE} \f$*/
	short has_ln_rcn_teee;/**< DLM: True if have \f$ C_l^{\phi\phi,TTTE} \f$*/
	short has_ln_rcn_tbeb;/**< DLM: True if have \f$ C_l^{\phi\phi,TTTE} \f$*/

	int nlt_size;  /**< number of noise \f$ N_l\f$ types requested */

	int index_nl_minvar; /**< DLM: index for type \f$ C_l^{\phi\phi,mv} \f$*/
	int index_nl_tt;  	 /**< DLM: index for type \f$ C_l^{\phi\phi,TT} \f$*/
	int index_nl_te; 	 /**< DLM: index for type \f$ C_l^{\phi\phi,TE} \f$*/
	int index_nl_ee;	 /**< DLM: index for type \f$ C_l^{\phi\phi,EE} \f$*/
	int index_nl_bb; 	 /**< DLM: index for type \f$ C_l^{\phi\phi,BB} \f$*/
	int index_nl_tb; 	 /**< DLM: index for type \f$ C_l^{\phi\phi,TB} \f$*/
	int index_nl_eb; 	 /**< DLM: index for type \f$ C_l^{\phi\phi,EB} \f$*/
	int index_nl_ttte; 	 /**< DLM: index for type \f$ C_l^{\phi\phi,TTTE} \f$*/
	int index_nl_ttee; 	 /**< DLM: index for type \f$ C_l^{\phi\phi,TTEE} \f$*/
	int index_nl_teee; 	 /**< DLM: index for type \f$ C_l^{\phi\phi,TEEE} \f$*/
	int index_nl_tbeb; 	 /**< DLM: index for type \f$ C_l^{\phi\phi,TBEB} \f$*/

	short calculate_derviaties_wrt_unlensed; /**< DLM:  */
	short lensed_wrt_unlensed; /**< DLM:  */
	short delensed_wrt_unlensed; /**< DLM:  */
    
	int l_tn_size, l_rcn_size, l_pn_size, l_cl_size, l_cl_lensed_size; /* DLM */

	double * l_cl; /**< DLM \f$*/
    double * l_cl_lensed; /**< DLM \f$*/

    short external_cl, external_cl_lensed; /* DLM */

    double * cl_tt_ext; /* DLM */
    double * cl_te_ext; /* DLM */
    double * cl_ee_ext; /* DLM */
    double * cl_bb_ext; /* DLM */
    double * cl_pp_ext; /* DLM */

    double * cl_lensed_tt_ext; /* DLM */
    double * cl_lensed_te_ext; /* DLM */
    double * cl_lensed_ee_ext; /* DLM */
    double * cl_lensed_bb_ext; /* DLM */
    double * cl_lensed_pp_ext; /* DLM */


	double sigma_beam; /**< DLM:  beamsize in radians. */
	double delta_noise; /**< DLM:  instrumental noise in muK-radians. */

	//@}

  /** @name - technical parameters */

  //@{

  short lensing_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */

	short delensing_verbose; /**< flag regulating the amount of information sent to standard output (none if set to zero) */ /* DLM */

  ErrorMsg error_message; /**< zone for writing error messages */

  //@}
};

/*************************************************************************************************************/
/* @cond INCLUDE_WITH_DOXYGEN */
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

	int delensing_cl_at_l(
						  struct lensing * ple,
						  int l,
						  double * cl_delensed
						  );

	int lensing_init(
					 struct precision * ppr,
					 struct perturbations * ppt,
					 struct harmonic * phr,
					 struct fourier * pfo,
					 struct lensing * ple
					 );

	int lensing_free(
					 struct lensing * ple
					 );

	int lensing_indices(
						struct precision * ppr,
						struct harmonic * phr,
						struct lensing * ple
						);

	int lensing_lensed_cl_tt(
							 double *ksi,
							 double **d00,
							 double *w8,
							 int nmu,
							 struct lensing * ple
							 );
    int lensing_lensed_cl_tt_derv(
                             double **ksi,
                             double **d00,
                             double *w8,
                             int nmu,
                             struct lensing * ple
                             );
    int lensing_lensed_cl_tt_derv_all(
                             double **ksi,
                             double **d00,
                             double *w8,
                             int nmu,
                             struct lensing * ple
                             );
    int lensing_delensed_cl_tt_derv(
                             double **ksi_dlu_derv,
                             double **d00,
                             double *w8,
                             int nmu,
                             struct lensing * ple
                             );
    int lensing_delensed_cl_tt_derv_all(
                             double **ksi_dlu_derv,
                             double **d00,
                             double *w8,
                             int nmu,
                             struct lensing * ple
                             );

	int lensing_delensed_cl_tt(
							   double *ksi_dl,
							   double *hla,
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
    int lensing_lensed_cl_te_derv(
                             double **ksiX,
                             double **d20,
                             double *w8,
                             int nmu,
                             struct lensing * ple
                             );
    int lensing_lensed_cl_te_derv_all(
                             double **ksiX,
                             double **d20,
                             double *w8,
                             int nmu,
                             struct lensing * ple
                             );
    int lensing_delensed_cl_te_derv(
                             double **ksiX_dlu_derv,
                             double **d20,
                             double *w8,
                             int nmu,
                             struct lensing * ple
                             );
    int lensing_delensed_cl_te_derv_all(
                             double **ksiX_dlu_derv,
                             double **d20,
                             double *w8,
                             int nmu,
                             struct lensing * ple
                             );
                             
	int lensing_delensed_cl_te(
							   double *ksiX_dl,
							   double *hla,
							   double *hla_P,
							   double **d20,
							   double *w8,
							   int nmu,
							   struct lensing * ple
							   );

	int lensing_lensed_cl_ee_bb(
								double *ksip,
								double *ksim,
								double **d22,
								double **d2m2,
								double *w8,
								int nmu,
								struct lensing * ple
								);

    int lensing_lensed_cl_ee_bb_dervE(
                                double **ksip,
                                double **ksim,
                                double **d22,
                                double **d2m2,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                );
    int lensing_lensed_cl_ee_bb_dervE_all(
                                double **ksip,
                                double **ksim,
                                double **d22,
                                double **d2m2,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                );
    int lensing_delensed_cl_ee_bb_dervE(
                                double **ksip_dlu_dervE,
                                double **ksim_dlu_dervE,
                                double **d22,
                                double **d2m2,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                );
    int lensing_delensed_cl_ee_bb_dervE_all(
                                double **ksip_dlu_dervE,
                                double **ksim_dlu_dervE,
                                double **d22,
                                double **d2m2,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                );

    int lensing_lensed_cl_ee_bb_dervB(
                                double **ksip,
                                double **ksim,
                                double **d22,
                                double **d2m2,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                );
    int lensing_lensed_cl_ee_bb_dervB_all(
                                double **ksip,
                                double **ksim,
                                double **d22,
                                double **d2m2,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                );
    int lensing_delensed_cl_ee_bb_dervB(
                                double **ksip_dlu_dervB,
                                double **ksim_dlu_dervB,
                                double **d22,
                                double **d2m2,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                );
    int lensing_delensed_cl_ee_bb_dervB_all(
                                double **ksip_dlu_dervB,
                                double **ksim_dlu_dervB,
                                double **d22,
                                double **d2m2,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                );

	int lensing_delensed_cl_ee_bb(
								  double *ksip_dl,
								  double *ksim_dl,
								  double *hla,
								  double **d22,
								  double **d2m2,
								  double *w8,
								  int nmu,
								  struct lensing * ple
								  );

	int lensing_addback_cl_tt(
							  struct lensing *ple,
							  double *cl_tt
							  );
//    int lensing_addback_cl_tt_derv(
//                              struct lensing *ple,
//                              double *cl_tt
//                              );

	int lensing_addback_cl_te(
							  struct lensing *ple,
							  double *cl_te
							  );
//    int lensing_addback_cl_te_derv(
//                              struct lensing *ple,
//                              double *cl_te
//                              );

	int lensing_addback_cl_ee_bb(
								 struct lensing *ple,
								 double *cl_ee,
								 double *cl_bb
								 );
//    int lensing_addback_cl_ee_bb_dervE(
//                                 struct lensing *ple,
//                                 double *cl_ee,
//                                 double *cl_bb
//                                 );
//    int lensing_addback_cl_ee_bb_dervB(
//                                 struct lensing *ple,
//                                 double *cl_ee,
//                                 double *cl_bb
//                                 );


	int lensing_addback_cl_dl_tt(
							  struct lensing *ple,
							  double *cl_tt
							  );

	int lensing_addback_cl_dl_te(
							  struct lensing *ple,
							  double *cl_te
							  );

	int lensing_addback_cl_dl_ee_bb(
								 struct lensing *ple,
								 double *cl_ee,
								 double *cl_bb
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

	int lensing_d01(
					double * mu,
					int num_mu,
					int lmax,
					double ** d01
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

	int lensing_d2m1(
					 double * mu,
					 int num_mu,
					 int lmax,
					 double ** d2m1
					 );

	int lensing_d21(
					double * mu,
					int num_mu,
					int lmax,
					double ** d21
					);

	int lensing_d20(
					double * mu,
					int num_mu,
					int lmax,
					double ** d20
					);

	int lensing_d30(
					double * mu,
					int num_mu,
					int lmax,
					double ** d30
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

	int lensing_d3m2(
					 double * mu,
					 int num_mu,
					 int lmax,
					 double ** d3m2
					 );

	int lensing_d32(
					double * mu,
					int num_mu,
					int lmax,
					double ** d32
					);

	int lensing_d3m3(
					 double * mu,
					 int num_mu,
					 int lmax,
					 double ** d3m3
					 );

	int lensing_d33(
					double * mu,
					int num_mu,
					int lmax,
					double ** d33
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

	int idealized_temperature_noise_spectrum_init(
												  struct lensing *ple,
												  struct harmonic * phr
												  );

	int idealized_polarization_noise_spectrum_init(
												   struct lensing *ple,
												   struct harmonic * phr
												   );

	int external_lens_recon_noise_spectrum_init(
												struct lensing *ple
												);

	int external_temperature_noise_spectrum_init(
												 struct lensing *ple,
												 struct harmonic * phr
												 );

	int external_polarization_noise_spectrum_init(
												  struct lensing *ple,
												  struct harmonic * phr
												  );
	int external_cmb_spectra_init(
												  struct lensing *ple,
												  struct harmonic * phr
												  );
	int external_lensed_cmb_spectra_init(
												  struct lensing *ple,
												  struct harmonic * phr
												  );

	int lensing_reconstruction_noise_tt(
										double **d00,
										double **d01,
										double **d11,
										double **d1m1,
										double *w8,
										int nmu,
										double *cl_tt_unlensed,
										double *cl_tt_lensed,
										struct lensing * ple,
										struct harmonic * phr);

	int lensing_reconstruction_noise_eb(
										double **d3m3,
										double **d33,
										double **d3m1,
										double **d31,
										double **d1m1,
										double **d11,
										double **d2m2,
										double **d22,
										double **d2m1,
										double **d21,
										double **d3m2,
										double **d32,
										double *w8,
										int nmu,
										double *cl_ee_th,
										double *cl_bb_th,
										double *cl_ee_ln,
										double *cl_bb_ln,
										struct lensing * ple,
										struct harmonic * phr
										);

	int lensing_reconstruction_noise_tb(
										double **d11,
										double **d1m1,
										double **d22,
										double **d2m2,
										double **d31,
										double **d3m1,
										double **d33,
										double **d3m3,
										double *w8,
										int nmu,
										double *cl_te_th,
										double *cl_tt_ln,
										double *cl_bb_ln,
										struct lensing * ple,
										struct harmonic * phr);

	int lensing_reconstruction_noise_te(
										double **d00,
										double **d01,
										double **d1m1,
										double **d11,
										double **d2m1,
										double **d21,
										double **d2m2,
										double **d22,
										double **d30,
										double **d3m1,
										double **d31,
										double **d3m3,
										double **d33,
										double *w8,
										int nmu,
										double *cl_te_th,
										double *cl_tt_ln,
										double *cl_ee_ln,
										struct lensing * ple,
										struct harmonic * phr
										);

	int lensing_reconstruction_noise_pp(
										double **d3m3,
										double **d33,
										double **d3m1,
										double **d31,
										double **d1m1,
										double **d11,
										double **d2m2,
										double **d22,
										double **d2m1,
										double **d21,
										double **d3m2,
										double **d32,
										double *w8,
										int nmu,
										double *cl_pol_th,
										double *cl_pol_ln,
										int ee_or_bb,
										struct lensing * ple,
										struct harmonic * phr
										);
	int lensing_reconst_nl_at_l(
								struct lensing * ple,
								int l,
								double * nl_rec
								);

	int lensing_reconstr_nl_tt(
							   double *zeta_t00,
							   double *zeta_t01,
							   double *zeta_t0m1,
							   double *zeta_t11,
							   double *zeta_t1m1,
							   double **d11,
							   double **d1m1,
							   double *w8,
							   int nmu,
							   struct lensing *ple);

	int lensing_reconstr_nl_minvar(
							   double **d3m3,
							   double **d33,
							   double **d3m2,
							   double **d32,
							   double **d3m1,
							   double **d31,
							   double **d30,
							   double **d2m2,
							   double **d22,
							   double **d2m1,
							   double **d21,
							   double **d20,
							   double **d1m1,
							   double **d11,
							   double **d01,
							   double **d00,
							   double *w8,
							   int nmu,
							   double *cl_tt_th,double *cl_tt_ln,
							   double *cl_te_th,double *cl_te_ln,
							   double *cl_ee_th,double *cl_ee_ln,
							   double *cl_bb_th,double *cl_bb_ln,
							   struct lensing * ple,
							   struct harmonic * phr
							   );

	int lensing_reconstr_nl_EB(
							   double **d3m3,  double **d33,
							   double **d3m2,  double **d32,
							   double **d3m1,  double **d31,
							   double **d2m2,  double **d22,
							   double **d2m1,  double **d21,
							   double **d1m1,  double **d11,
							   double **d00,
							   double *w8,
							   int nmu,
							   double *cl_ee_th,double *cl_ee_ln,
							   double *cl_bb_th,double *cl_bb_ln,
							   struct lensing * ple,
							   struct harmonic * phr
							   );

	int dl_spectra_bb_ala_smith(
								double **d3m3,  double **d33,
								double **d3m2,  double **d32,
								double **d3m1,  double **d31,
								double **d2m2,  double **d22,
								double **d1m1,  double **d11,
								double **d00,
								double *w8,
								int nmu,
								double *cl_ee_th,double *cl_ee_ln,
								double *cl_pp,double *nl_ebeb,
								struct lensing * ple,
								struct harmonic * phr
								);


	int lensing_delensed_derv_cl_tt(
									double **ksi_dl_derv,
									double *hla,
									double **d00,
									double *w8,
									int nmu,
									struct lensing * ple
									);

	int lensing_addback_derv_cl_dl_tt(
									  struct lensing * ple,
									  double *cl_tt);

	int lensing_delensed_derv_cl_te(
						   double **ksiX_dl_derv,
						   double *hla,
						   double *hla_P,
						   double **d20,
						   double *w8,
						   int nmu,
						   struct lensing * ple
						   );

	int lensing_addback_derv_cl_dl_te(
									  struct lensing * ple,
									  double *cl_te);

	int lensing_delensed_derv_cl_ee_bb(
									   double **ksip_dl_derv,
									   double **ksim_dl_derv,
									   double *hla,
									   double **d22,
									   double **d2m2,
									   double *w8,
									   int nmu,
									   struct lensing * ple
									   );

	int lensing_addback_derv_cl_dl_ee_bb(
										 struct lensing * ple,
										 double * cl_ee,
										 double * cl_bb);


	int delensing_derv_cl_at_l(
							   struct lensing * ple,
							   int l,
							   int k,
							   double * cl_delensed
							   );

	int delensing_derv_cl_tt_at_l(
								  struct lensing * ple,
								  double k,
								  double l,
								  double * cl_derv_delensed);

    int lensing_derv_cl_tt_tt_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_lens_derv_TT_TT);

    int lensing_derv_dl_cl_tt_tt_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_TT_TT);

	int delensing_derv_cl_te_at_l( /* DLM */
								  struct lensing * ple,
								  double k,
								  double l,
								  double * cl_derv_delensed);

    int lensing_derv_cl_te_te_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_lens_derv_TE_TE);

    int lensing_derv_dl_cl_te_te_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_TE_TE);

	int delensing_derv_cl_ee_at_l( /* DLM */
								  struct lensing * ple,
								  double k,
								  double l,
								  double * cl_derv_delensed);

    int lensing_derv_cl_ee_ee_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_lens_derv_EE_EE);

    int lensing_derv_dl_cl_ee_ee_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_EE_EE);

    int lensing_derv_cl_ee_bb_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_lens_derv_EE_BB);

    int lensing_derv_dl_cl_ee_bb_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_EE_BB);

    int lensing_derv_cl_bb_ee_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_lens_derv_BB_EE);

    int lensing_derv_dl_cl_bb_ee_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_BB_EE);

    int lensing_derv_cl_bb_bb_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_lens_derv_BB_BB);

    int lensing_derv_dl_cl_bb_bb_at_l(
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_BB_BB);

	int delensing_derv_cl_bb_at_l( /* DLM */
								  struct lensing * ple,
								  double k,
								  double l,
								  double * cl_derv_delensed);

#ifdef __cplusplus
}
#endif

#endif
/* @endcond */
