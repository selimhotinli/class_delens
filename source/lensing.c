/** @file lensing.c Documented lensing module
 *
 * Simon Prunet and Julien Lesgourgues, 6.12.2010
 *
 * This module computes the lensed temperature and polarization
 * anisotropy power spectra \f$ C_l^{X}, P(k), ... \f$'s given the
 * unlensed temperature, polarization and lensing potential spectra.
 *
 * Follows Challinor and Lewis full-sky method, astro-ph/0502425
 *
 * DLM: Follows Green, Meyers and van Engelen, arvix:1609.08143
 *
 * The following functions can be called from other modules:
 *
 * -# lensing_init() at the beginning (but after harmonic_init())
 * -# lensing_cl_at_l() at any time for computing Cl_lensed at any l
 * -# lensing_free() at the end
 */

#include "lensing.h"
#include <time.h>
#include <math.h> /* DLM */

/**
 * Anisotropy power spectra \f$ C_l\f$'s for all types, modes and initial conditions.
 * SO FAR: ONLY SCALAR
 *
 * This routine evaluates all the lensed \f$ C_l\f$'s at a given value of l by
 * picking it in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 *
 * This function can be called from whatever module at whatever time,
 * provided that lensing_init() has been called before, and
 * lensing_free() has not been called yet.
 *
 * @param ple        Input: pointer to lensing structure
 * @param l          Input: multipole number
 * @param cl_lensed  Output: lensed \f$ C_l\f$'s for all types (TT, TE, EE, etc..)
 * @return the error status
 */

int lensing_cl_at_l(
                    struct lensing * ple,
                    int l,
                    double * cl_lensed    /* array with argument cl_lensed[index_ct] (must be already allocated) */
) {
    int last_index;
    int index_lt;
    
    class_test(l > ple->l_lensed_max,
               ple->error_message,
               "you asked for lensed Cls at l=%d, they were computed only up to l=%d, you should increase l_max_scalars or decrease the precision parameter delta_l_max",l,ple->l_lensed_max);
    
    class_call(array_interpolate_spline(ple->l,
                                        ple->l_size,
                                        ple->cl_lens,
                                        ple->ddcl_lens,
                                        ple->lt_size,
                                        l,
                                        &last_index,
                                        cl_lensed,
                                        ple->lt_size,
                                        ple->error_message),
               ple->error_message,
               ple->error_message);
    
    /* set to zero for the types such that l<l_max */
    for (index_lt=0; index_lt<ple->lt_size; index_lt++)
    if ((int)l > ple->l_max_lt[index_lt])
        cl_lensed[index_lt]=0.;
    
    return _SUCCESS_;
}

int delensing_cl_at_l( /* DLM */
                      struct lensing * ple,
                      int l,
                      double * cl_delensed    /* array with argument cl_delens[index_ct] (must be already allocated) */
) {
    int last_index;
    int index_lt;
    
    class_test(l > ple->l_lensed_max,
               ple->error_message,
               "you asked for lensed Cls at l=%d, they were computed only up to l=%d, you should increase l_max_scalars or decrease the precision parameter delta_l_max",l,ple->l_lensed_max);
    
    class_call(array_interpolate_spline(ple->l_dl,
                                        ple->dl_size,
                                        ple->cl_delens,
                                        ple->ddcl_delens,
                                        ple->dlt_size,
                                        l,
                                        &last_index,
                                        cl_delensed,
                                        ple->dlt_size,
                                        ple->error_message),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}

// derivative of delensed TT w.r.t. lensing potential
int delensing_derv_cl_tt_at_l( /* DLM */
                              struct lensing * ple,
                              double k,
                              double l,
                              double * cl_derv_delensed) {
    int last_index;
    int index_lt;
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    
    class_call(dlm_splin2(ple->l_dl,
                          ple->l_dl,
                          ple->cl_dl_tt_pderv,
                          ple->ddcl_dl_tt_pderv,
                          ple->dl_size-1,
                          ple->dl_size-1,
                          k, // L = Cl^\phi\phi wavenumber
                          l, // l = Cl^TT wavenumber
                          cl_derv_delensed),
               ple->error_message,
               ple->error_message);
    return _SUCCESS_;
}
// derivative of lensed TT w.r.t. unlensed TT
int lensing_derv_cl_tt_tt_at_l( /* DLM */
                               struct lensing * ple,
                               double k,
                               double l,
                               double * cl_lens_derv_TT_TT) {
    int last_index;
    int index_lt;
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_lens_derv_TT_TT,
                          ple->ddcl_lens_derv_TT_TT,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_lens_derv_TT_TT),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of delensed TT w.r.t. unlensed TT
int lensing_derv_dl_cl_tt_tt_at_l( /* DLM */
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_TT_TT) {
    int last_index;
    int index_lt;
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_delens_derv_TT_TT,
                          ple->ddcl_delens_derv_TT_TT,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for delensed CL
                          l, // l for unlensed CL
                          cl_delens_derv_TT_TT),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of delensed TE w.r.t. lensing potential
int delensing_derv_cl_te_at_l( /* DLM */
                              struct lensing * ple,
                              double k,
                              double l,
                              double * cl_derv_delensed) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l_dl,
                          ple->l_dl,
                          ple->cl_dl_te_pderv,
                          ple->ddcl_dl_te_pderv,
                          ple->dl_size-1,
                          ple->dl_size-1,
                          k, // L = Cl^\phi\phi wavenumber
                          l, // l = Cl^TT wavenumber
                          cl_derv_delensed),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of lensed TE w.r.t. unlensed TE
int lensing_derv_cl_te_te_at_l( /* DLM */
                               struct lensing * ple,
                               double k,
                               double l,
                               double * cl_lens_derv_TE_TE) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_lens_derv_TE_TE,
                          ple->ddcl_lens_derv_TE_TE,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_lens_derv_TE_TE),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of delensed TE w.r.t. unlensed TE
int lensing_derv_dl_cl_te_te_at_l( /* DLM */
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_TE_TE) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_delens_derv_TE_TE,
                          ple->ddcl_delens_derv_TE_TE,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_delens_derv_TE_TE),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of delensed EE w.r.t. lensing potential
int delensing_derv_cl_ee_at_l( /* DLM */
                              struct lensing * ple,
                              double k,
                              double l,
                              double * cl_derv_delensed) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l_dl,
                          ple->l_dl,
                          ple->cl_dl_ee_pderv,
                          ple->ddcl_dl_ee_pderv,
                          ple->dl_size-1,
                          ple->dl_size-1,
                          k, // L = Cl^\phi\phi wavenumber
                          l, // l = Cl^TT wavenumber
                          cl_derv_delensed),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}

// derivative of lensed EE w.r.t. unlensed EE
int lensing_derv_cl_ee_ee_at_l( /* DLM */
                               struct lensing * ple,
                               double k,
                               double l,
                               double * cl_lens_derv_EE_EE) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_lens_derv_EE_EE,
                          ple->ddcl_lens_derv_EE_EE,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_lens_derv_EE_EE),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of delensed EE w.r.t. unlensed EE
int lensing_derv_dl_cl_ee_ee_at_l( /* DLM */
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_EE_EE) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_delens_derv_EE_EE,
                          ple->ddcl_delens_derv_EE_EE,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_delens_derv_EE_EE),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}

// derivative of lensed EE w.r.t. unlensed BB
int lensing_derv_cl_ee_bb_at_l( /* DLM */
                               struct lensing * ple,
                               double k,
                               double l,
                               double * cl_lens_derv_EE_BB) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_lens_derv_EE_BB,
                          ple->ddcl_lens_derv_EE_BB,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_lens_derv_EE_BB),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of delensed EE w.r.t. unlensed BB
int lensing_derv_dl_cl_ee_bb_at_l( /* DLM */
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_EE_BB) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_delens_derv_EE_BB,
                          ple->ddcl_delens_derv_EE_BB,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_delens_derv_EE_BB),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}

// derivative of lensed BB w.r.t. unlensed EE
int lensing_derv_cl_bb_ee_at_l( /* DLM */
                               struct lensing * ple,
                               double k,
                               double l,
                               double * cl_lens_derv_BB_EE) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_lens_derv_BB_EE,
                          ple->ddcl_lens_derv_BB_EE,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_lens_derv_BB_EE),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of delensed BB w.r.t. unlensed EE
int lensing_derv_dl_cl_bb_ee_at_l( /* DLM */
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_BB_EE) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_delens_derv_BB_EE,
                          ple->ddcl_delens_derv_BB_EE,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_delens_derv_BB_EE),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}
// derivative of lensed BB w.r.t. unlensed BB
int lensing_derv_cl_bb_bb_at_l( /* DLM */
                               struct lensing * ple,
                               double k,
                               double l,
                               double * cl_lens_derv_BB_BB) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_lens_derv_BB_BB,
                          ple->ddcl_lens_derv_BB_BB,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_lens_derv_BB_BB),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}

// derivative of delensed BB w.r.t. unlensed BB
int lensing_derv_dl_cl_bb_bb_at_l( /* DLM */
                                  struct lensing * ple,
                                  double k,
                                  double l,
                                  double * cl_delens_derv_BB_BB) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l,
                          ple->l,
                          ple->cl_delens_derv_BB_BB,
                          ple->ddcl_delens_derv_BB_BB,
                          ple->l_size-1,
                          ple->l_size-1,
                          k, // l for lensed CL
                          l, // l for unlensed CL
                          cl_delens_derv_BB_BB),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}

// derivative of delensed BB w.r.t. lensing potential
int delensing_derv_cl_bb_at_l( /* DLM */
                              struct lensing * ple,
                              double k,
                              double l,
                              double * cl_derv_delensed) {
    int last_index;
    int index_lt;
    
    /*This routine returns an interpolated function value cl_derv_delensed
     by bicubic spline interpolation at k and l.*/
    class_call(dlm_splin2(ple->l_dl,
                          ple->l_dl,
                          ple->cl_dl_bb_pderv,
                          ple->ddcl_dl_bb_pderv,
                          ple->dl_size-1,
                          ple->dl_size-1,
                          k, // L = Cl^\phi\phi wavenumber
                          l, // l = Cl^TT wavenumber
                          cl_derv_delensed),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}


int lensing_reconst_nl_at_l( /* DLM */
                            struct lensing * ple,
                            int l,
                            double * nl_rec    /* array with argument cl_delens[index_ct] (must be already allocated) */
) {
    int last_index;
    int index_lt;
    
    class_test(l > ple->l_lensed_max,
               ple->error_message,
               "you asked for lensed Cls at l=%d, they were computed only up to l=%d, you should increase l_max_scalars or decrease the precision parameter delta_l_max",l,ple->l_lensed_max);
    
    class_call(array_interpolate_spline(ple->l,
                                        ple->dl_size,
                                        ple->nl_rcn,
                                        ple->ddnl_rcn,
                                        ple->nlt_size,
                                        l,
                                        &last_index,
                                        nl_rec,
                                        ple->nlt_size,
                                        ple->error_message),
               ple->error_message,
               ple->error_message);
    
    return _SUCCESS_;
}



/**
 * This routine initializes the lensing structure (in particular,
 * computes table of lensed anisotropy spectra \f$ C_l^{X} \f$)
 *
 * @param ppr Input: pointer to precision structure
 * @param ppt Input: pointer to perturbation structure (just in case, not used in current version...)
 * @param phr Input: pointer to harmonic structure
 * @param pfo Input: pointer to fourier structure
 * @param ple Output: pointer to initialized lensing structure
 * @return the error status
 */

int lensing_init(
                 struct precision * ppr,
                 struct perturbations * ppt,
                 struct harmonic * phr,
                 struct fourier * pfo,
                 struct lensing * ple
                 ) {
    
    /** Summary: */
    /** - Define local variables */
    
    double * mu; /* mu[index_mu]: discretized values of mu
                  between -1 and 1, roots of Legendre polynomial */
    double * w8; /* Corresponding Gauss-Legendre quadrature weights */
    double theta,delta_theta;
    
    double ** d00;  /* dmn[index_mu][index_l] */
    double ** d11;
    double ** d2m2;
    double ** d22 = NULL;
    double ** d20 = NULL;
    double ** d1m1;
    double ** d31 = NULL;
    double ** d40 = NULL;
    double ** d3m1 = NULL;
    double ** d3m3 = NULL;
    double ** d4m2 = NULL;
    double ** d4m4 = NULL;
    
    double * buf_dxx; /* buffer */
    double * buf2_dxx; /* buffer2 */
    double * buf3_dxx; /* buffer3 */
    double * buf4_dxx; /* buffer3 */
    
    double ** d33 = NULL; /* DLM */
    double ** d30 = NULL; /* DLM */
    double ** d01 = NULL; /* DLM */
    double ** d21 = NULL; /* DLM */
    double ** d2m1 = NULL; /* DLM */
    double ** d32  = NULL; /* DLM */
    double ** d3m2 = NULL; /* DLM */
    
    double * Cgl;   /* Cgl[index_mu] */
    
    double * Cgl_obs = NULL; double * Cgl_cross = NULL; /* DLM */
    
    double * Cgl2;  /* Cgl2[index_mu] */
    
    double * Cgl2_obs = NULL; double * Cgl2_cross = NULL; /* DLM */
    
    double * sigma2; /* sigma[index_mu] */
    
    double * sigma2_hhbar = NULL; double * sigma2_hh = NULL; /* DLM */
    
    double * ksi = NULL;   /* ksi[index_mu] */
    double * ksiX = NULL;  /* ksiX[index_mu] */
    double * ksip = NULL;  /* ksip[index_mu] */
    double * ksim = NULL;  /* ksim[index_mu] */
    
    double ** ksi_ln_derv = NULL;  /* ksi_ln_derv[index_mu] */
    double ** ksiX_ln_derv = NULL; /* ksiX_ln_derv[index_mu] */
    double ** ksip_ln_dervE = NULL; /* ksip_ln_dervE[index_mu] */
    double ** ksim_ln_dervE = NULL; /* ksim_ln_dervE[index_mu] */
    double ** ksip_ln_dervB = NULL; /* ksip_ln_dervB[index_mu] */
    double ** ksim_ln_dervB = NULL; /* ksim_ln_dervB[index_mu] */
    
    double ** ksi_dlu_derv = NULL;  /* ksi_dl_derv[index_mu] */
    double ** ksiX_dlu_derv = NULL; /* ksiX_dl_derv[index_mu] */
    double ** ksip_dlu_dervE = NULL; /* ksip_dl_dervE[index_mu] */
    double ** ksim_dlu_dervE = NULL; /* ksim_dl_dervE[index_mu] */
    double ** ksip_dlu_dervB = NULL; /* ksip_dl_dervB[index_mu] */
    double ** ksim_dlu_dervB = NULL; /* ksim_dl_dervB[index_mu] */
    
    double * ksi_dl = NULL;  double * ksiX_dl = NULL; /* DLM */
    double * ksip_dl = NULL; double * ksim_dl = NULL; /* DLM */
    
    double fac,fac1;
    
    double X_000;
    
    double X_000_hhbar; double X_000_hh; /* DLM */
    
    double X_p000;
    
    double X_p000_hhbar; double X_p000_hh; /* DLM */
    
    double X_220;
    
    double X_220_hhbar; double X_220_hh; /* DLM */
    
    double X_022;
    double X_p022;
    double X_121;
    double X_132;
    double X_242;
    
    double X_022_hhbar; double X_022_hh; /* DLM */
    double X_121_hhbar; double X_121_hh; /* DLM */
    double X_132_hhbar; double X_132_hh; /* DLM */
    double X_242_hhbar; double X_242_hh; /* DLM */
    
    double X_p022_hhbar;
    double X_p022_hh;
    
    /* comment what is below asap. */
    double facd;
    double sigma2_derv,sigma2_hh_derv,sigma2_hhbar_derv;
    double X_000_derv,X_000_hhbar_derv,X_000_hh_derv;
    double X_p000_derv,X_p000_hhbar_derv,X_p000_hh_derv;
    double X_022_derv,X_022_hhbar_derv,X_p022_hh_derv;
    double X_p022_derv,X_p022_hhbar_derv,X_022_hh_derv;
    double X_220_derv,X_220_hhbar_derv,X_220_hh_derv;
    double X_121_derv,X_121_hhbar_derv,X_121_hh_derv;
    double X_132_derv,X_132_hhbar_derv,X_132_hh_derv;
    double X_242_derv,X_242_hhbar_derv,X_242_hh_derv;
    
    double ** ksi_dl_derv;
    double ** ksiX_dl_derv;
    double ** ksip_dl_derv;
    double ** ksim_dl_derv;
    
    double    cl_delensed_derv;
    
    long int num_mu,index_mu,icount; /* DLM: 'long' int needed for \ell higher than ~12k when allocating d-arrays.*/
    int l,k,index_k,index_l,index1_l,index_ll,ii;
    double ll;
    
    double cratio; /* DLM */
    double mvold,tempd; /* DLM */
    
    double cratioEB; /* DLM */
    double mvoldEB,tempdEB; /* DLM */
    
    double * cratios; /* DLM */
    double * mvolds; double * tempds; /* DLM */
    
    double * cl_unlensed;  /* cl_unlensed[index_ct] */
    double * cl_tt; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
    double * cl_te = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
    double * cl_ee = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
    double * cl_bb = NULL; /* unlensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
    double * cl_pp; /* potential cl, to be filled to avoid repeated calls to harmonic_at_l */
    
    
    double * cl_pp_obs = NULL;    /* DLM: c_pp_obs = c_pp + nl_noise (observed lensing potential) */
    double * cl_pp_cross = NULL;  /* DLM: c_pp_cross = c_pp */
    
    double * gl = NULL; /* DLM: gl: delensing lensing potential filter  */
    double * hl = NULL; /* DLM: hl: delensing temperature filter 1/2 */
    double * hlbar = NULL;  /* DLM: hlbar: delensing temperature filter 2/2 */
    
    double * hl_P = NULL; /* DLM: hl_P: delensing polarization filter 1/2 */
    double * hlbar_P = NULL; /* DLM: hlbar_P: delensing polarization filter 2/2 */
    
    double res,resX,lens;
    double resp, resm, lensp, lensm;
    
    double res_dl_ml,resX_dl,resm_dl,resp_dl; /* DLM */
    double lens_dl,lensp_dl,lensm_dl; /* DLM */
    
    double res_dl_derv,resX_dl_derv,resm_dl_derv,resp_dl_derv; /* DLM */
    double lens_dl_derv,lensX_dl_derv,lensm_dl_derv,lensp_dl_derv; /* DLM */
    
    double * cl_lensed = NULL;
    
    double * cl_lens_tt = NULL;/* lensed  cl, to be filled to avoid repeated calls to lensing_cl_at_l */
    double * cl_lens_te = NULL;/* lensed  cl, to be filled to avoid repeated calls to lensing_cl_at_l */
    double * cl_lens_ee = NULL;/* lensed  cl, to be filled to avoid repeated calls to lensing_cl_at_l */
    double * cl_lens_bb = NULL;/* lensed  cl, to be filled to avoid repeated calls to lensing_cl_at_l */
    
    double * cl_delensed = NULL; /* DLM: similar purpose to cl_lensed */
    
    double * cl_delens_tt = NULL;/* DLM: delensed  cl, to be filled to avoid repeated calls to delensing_cl_at_l during iterative delensing */
    double * cl_delens_te = NULL;/* DLM: delensed  cl, to be filled to avoid repeated calls to delensing_cl_at_l during iterative delensing */
    double * cl_delens_ee = NULL;/* DLM: delensed  cl, to be filled to avoid repeated calls to delensing_cl_at_l during iterative delensing */
    double * cl_delens_bb = NULL;/* DLM: delensed  cl, to be filled to avoid repeated calls to delensing_cl_at_l during iterative delensing */
    
    double * nl_lensed = NULL; /* DLM: similar purpose to cl_lensed */
    
    double * min_varr_ln = NULL; /* DLM: minimum varriance lensing noise estimate */
    double * min_ebeb_ln = NULL; /* DLM: the EBEB varriance lensing noise estimate */
    
    int itr_index; /* DLM */
    
    double * sqrt1;
    double * sqrt2;
    double * sqrt3;
    double * sqrt4;
    double * sqrt5;
    
    double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */
    
    double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */
    
    double kk; // DLM.
    
    short type_2_flag = _FALSE_;
    short print_conv  = _FALSE_; // DLM: Flag for announcing that the minimum varriance spectrum has converged.
    
    int index_md;
    
    /* Timing */
    //double debut, fin;
    //double cpu_time;
    
    /** - check that we really want to compute at least one spectrum */
    
    if (ple->has_lensed_cls == _FALSE_) {
        if (ple->lensing_verbose > 0)
            printf("No lensing requested. Lensing module skipped.\n");
        return _SUCCESS_;
    }
    else {
        if (ple->lensing_verbose > 0) {
            printf("Computing lensed spectra ");
            if (ppr->accurate_lensing==_TRUE_)
                printf("(accurate mode)\n");
            else
                printf("(fast mode)\n");
        }
    }
    
    /** - initialize indices and allocate some of the arrays in the
     lensing structure */
    
    class_call(lensing_indices(ppr,phr,ple),
               ple->error_message,
               ple->error_message);
    
    /** - put all precision variables hare; will be stored later in precision structure */
    /** - Last element in \f$ \mu \f$ will be for \f$ \mu=1 \f$, needed for sigma2.
     The rest will be chosen as roots of a Gauss-Legendre quadrature **/
    
    if (ppr->accurate_lensing == _TRUE_) {
        num_mu=(ple->l_unlensed_max+ppr->num_mu_minus_lmax); /* Must be even ?? CHECK */
        num_mu += num_mu%2; /* Force it to be even */
    } else {
        /* Integrate correlation function difference on [0,pi/16] */
        num_mu = (ple->l_unlensed_max * 2 )/16;
    }
  }

  /** - Compute \f$ d^l_{mm'} (\mu) \f$*/

  icount = 0;
  class_alloc(d00,
              num_mu*sizeof(double*),
              ple->error_message);

  class_alloc(d11,
              num_mu*sizeof(double*),
              ple->error_message);

  class_alloc(d1m1,
              num_mu*sizeof(double*),
              ple->error_message);

  class_alloc(d2m2,
              num_mu*sizeof(double*),
              ple->error_message);
  icount += 4*num_mu*(ple->l_unlensed_max+1);

  if (ple->has_te==_TRUE_) {

    class_alloc(d20,
                num_mu*sizeof(double*),
                ple->error_message);
    /* Reserve last element of mu for mu=1, needed for sigma2 */
    mu[num_mu-1] = 1.0;
    
    class_alloc(w8,
                (num_mu-1)*sizeof(double),
                ple->error_message);
    
    if (ppr->accurate_lensing == _TRUE_) {
        
        //debut = omp_get_wtime();
        class_call(quadrature_gauss_legendre(mu,
                                             w8,
                                             num_mu-1,
                                             ppr->tol_gauss_legendre,
                                             ple->error_message),
                   ple->error_message,
                   ple->error_message);
        //fin = omp_get_wtime();
        //cpu_time = (fin-debut);
        //printf("time in quadrature_gauss_legendre=%4.3f s\n",cpu_time);
        
    } else { /* Crude integration on [0,pi/16]: Riemann sum on theta */
        
        delta_theta = _PI_/16. / (double)(num_mu-1);
        
        for (index_mu=0;index_mu<num_mu-1;index_mu++) {
            theta = (index_mu+1)*delta_theta;
            mu[index_mu] = cos(theta);
            w8[index_mu] = sin(theta)*delta_theta; /* We integrate on mu */
            
        }
    }
    
    /** - Compute \f$ d^l_{mm'} (\mu) \f$*/
    
    icount = 0;
    class_alloc(d00,
                num_mu*sizeof(double*),
                ple->error_message);
    
    class_alloc(d11,
                num_mu*sizeof(double*),
                ple->error_message);
    
    class_alloc(d1m1,
                num_mu*sizeof(double*),
                ple->error_message);
    
    class_alloc(d2m2,
                num_mu*sizeof(double*),
                ple->error_message);
    icount += 4*num_mu*(ple->l_unlensed_max+1);
    
    if(ple->has_te==_TRUE_) {
        
        class_alloc(d20,
                    num_mu*sizeof(double*),
                    ple->error_message);
        
        class_alloc(d3m1,
                    num_mu*sizeof(double*),
                    ple->error_message);
        
        class_alloc(d4m2,
                    num_mu*sizeof(double*),
                    ple->error_message);
        icount += 3*num_mu*(ple->l_unlensed_max+1);
    }
    
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        
        class_alloc(d22,
                    num_mu*sizeof(double*),
                    ple->error_message);
        
        class_alloc(d31,
                    num_mu*sizeof(double*),
                    ple->error_message);
        
        class_alloc(d3m3,
                    num_mu*sizeof(double*),
                    ple->error_message);
        
        if (ple->has_delensed_cls == _TRUE_){
            
            class_alloc(d33,
                        num_mu*sizeof(double*),
                        ple->error_message); /* DLM */
            
            class_alloc(d01,
                        num_mu*sizeof(double*),
                        ple->error_message); /* DLM */
            
            class_alloc(d21,
                        num_mu*sizeof(double*),
                        ple->error_message); /* DLM */
            
            class_alloc(d2m1,
                        num_mu*sizeof(double*),
                        ple->error_message); /* DLM */
            
            class_alloc(d32,
                        num_mu*sizeof(double*),
                        ple->error_message); /* DLM */
            
            class_alloc(d3m2,
                        num_mu*sizeof(double*),
                        ple->error_message); /* DLM */
            
            class_alloc(d30,
                        num_mu*sizeof(double*),
                        ple->error_message); /* DLM */
        }
        
        class_alloc(d40,
                    num_mu*sizeof(double*),
                    ple->error_message);
        
        class_alloc(d4m4,
                    num_mu*sizeof(double*),
                    ple->error_message);
        
        icount += 5*num_mu*(ple->l_unlensed_max+1);
        
        if (ple->has_delensed_cls == _TRUE_) icount += 7*num_mu*(ple->l_unlensed_max+1); /* DLM */
        
    }
    
    icount += 5*(ple->l_unlensed_max+1); /* for arrays sqrt1[l] to sqrt5[l] */
    
    /** - Allocate main contiguous buffer **/
    class_alloc(buf_dxx,
                icount * sizeof(double),
                ple->error_message);
    
    icount = 0;
    for (index_mu=0; index_mu<num_mu; index_mu++) {
        
        d00[index_mu] = &(buf_dxx[icount+index_mu            * (ple->l_unlensed_max+1)]);
        d11[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)   * (ple->l_unlensed_max+1)]);
        d1m1[index_mu]= &(buf_dxx[icount+(index_mu+2*num_mu) * (ple->l_unlensed_max+1)]);
        d2m2[index_mu]= &(buf_dxx[icount+(index_mu+3*num_mu) * (ple->l_unlensed_max+1)]);
    }
    icount += 4*num_mu*(ple->l_unlensed_max+1);
    
    if (ple->has_te==_TRUE_) {
        for (index_mu=0; index_mu<num_mu; index_mu++) {
            d20[index_mu] = &(buf_dxx[icount+index_mu            * (ple->l_unlensed_max+1)]);
            d3m1[index_mu]= &(buf_dxx[icount+(index_mu+num_mu)   * (ple->l_unlensed_max+1)]);
            d4m2[index_mu]= &(buf_dxx[icount+(index_mu+2*num_mu) * (ple->l_unlensed_max+1)]);
        }
        icount += 3*num_mu*(ple->l_unlensed_max+1);
    }
    
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        
        for (index_mu=0; index_mu<num_mu; index_mu++) {
            d22[index_mu] = &(buf_dxx[icount+index_mu            * (ple->l_unlensed_max+1)]);
            d31[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)   * (ple->l_unlensed_max+1)]);
            d3m3[index_mu]= &(buf_dxx[icount+(index_mu+2*num_mu) * (ple->l_unlensed_max+1)]);
            
            if (ple->has_delensed_cls == _TRUE_){
                
                d33[index_mu]  = &(buf_dxx[icount+(index_mu+3*num_mu) * (ple->l_unlensed_max+1)]); /* DLM */
                d01[index_mu]  = &(buf_dxx[icount+(index_mu+4*num_mu) * (ple->l_unlensed_max+1)]); /* DLM */
                d21[index_mu]  = &(buf_dxx[icount+(index_mu+5*num_mu) * (ple->l_unlensed_max+1)]); /* DLM */
                d2m1[index_mu] = &(buf_dxx[icount+(index_mu+6*num_mu) * (ple->l_unlensed_max+1)]); /* DLM */
                d32[index_mu]  = &(buf_dxx[icount+(index_mu+7*num_mu) * (ple->l_unlensed_max+1)]); /* DLM */
                d3m2[index_mu] = &(buf_dxx[icount+(index_mu+8*num_mu) * (ple->l_unlensed_max+1)]); /* DLM */
                d30[index_mu]  = &(buf_dxx[icount+(index_mu+9*num_mu) * (ple->l_unlensed_max+1)]); /* DLM */
                d40[index_mu]  = &(buf_dxx[icount+(index_mu+10*num_mu) * (ple->l_unlensed_max+1)]);
                d4m4[index_mu] = &(buf_dxx[icount+(index_mu+11*num_mu) * (ple->l_unlensed_max+1)]);
                
            }
            
            else
            {
                d40[index_mu]  = &(buf_dxx[icount+(index_mu+3*num_mu) * (ple->l_unlensed_max+1)]);
                d4m4[index_mu] = &(buf_dxx[icount+(index_mu+4*num_mu) * (ple->l_unlensed_max+1)]);
            }
            
        }
        icount += 5*num_mu*(ple->l_unlensed_max+1);
        
        if (ple->has_delensed_cls == _TRUE_) icount += 7*num_mu*(ple->l_unlensed_max+1);
    }
    
    sqrt1 = &(buf_dxx[icount]);
    icount += ple->l_unlensed_max+1;
    sqrt2 = &(buf_dxx[icount]);
    icount += ple->l_unlensed_max+1;
    sqrt3 = &(buf_dxx[icount]);
    icount += ple->l_unlensed_max+1;
    sqrt4 = &(buf_dxx[icount]);
    icount += ple->l_unlensed_max+1;
    sqrt5 = &(buf_dxx[icount]);
    icount += ple->l_unlensed_max+1;
    
    //debut = omp_get_wtime();
    class_call(lensing_d00(mu,num_mu,ple->l_unlensed_max,d00),
               ple->error_message,
               ple->error_message);
    
    class_call(lensing_d11(mu,num_mu,ple->l_unlensed_max,d11),
               ple->error_message,
               ple->error_message);
    
    class_call(lensing_d1m1(mu,num_mu,ple->l_unlensed_max,d1m1),
               ple->error_message,
               ple->error_message);
    
    class_call(lensing_d2m2(mu,num_mu,ple->l_unlensed_max,d2m2),
               ple->error_message,
               ple->error_message);
    
    //fin = omp_get_wtime();
    //cpu_time = (fin-debut);
    //printf("time in lensing_dxx=%4.3f s\n",cpu_time);
    
    if (ple->has_te==_TRUE_) {
        
        class_call(lensing_d20(mu,num_mu,ple->l_unlensed_max,d20),
                   ple->error_message,
                   ple->error_message);
        
        class_call(lensing_d3m1(mu,num_mu,ple->l_unlensed_max,d3m1),
                   ple->error_message,
                   ple->error_message);
        
        class_call(lensing_d4m2(mu,num_mu,ple->l_unlensed_max,d4m2),
                   ple->error_message,
                   ple->error_message);
        
    }
    
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        
        class_call(lensing_d22(mu,num_mu,ple->l_unlensed_max,d22),
                   ple->error_message,
                   ple->error_message);
        
        class_call(lensing_d31(mu,num_mu,ple->l_unlensed_max,d31),
                   ple->error_message,
                   ple->error_message);
        
        class_call(lensing_d3m3(mu,num_mu,ple->l_unlensed_max,d3m3),
                   ple->error_message,
                   ple->error_message);
        
        if (ple->has_delensed_cls == _TRUE_){
            
            class_call(lensing_d33(mu,num_mu,ple->l_unlensed_max,d33),
                       ple->error_message,
                       ple->error_message); /* DLM */
            
            class_call(lensing_d01(mu,num_mu,ple->l_unlensed_max,d01),
                       ple->error_message,
                       ple->error_message); /* DLM */
            
            class_call(lensing_d21(mu,num_mu,ple->l_unlensed_max,d21),
                       ple->error_message,
                       ple->error_message); /* DLM */
            
            class_call(lensing_d2m1(mu,num_mu,ple->l_unlensed_max,d2m1),
                       ple->error_message,
                       ple->error_message); /* DLM */
            
            class_call(lensing_d30(mu,num_mu,ple->l_unlensed_max,d30),
                       ple->error_message,
                       ple->error_message); /* DLM */
            
            class_call(lensing_d32(mu,num_mu,ple->l_unlensed_max,d32),
                       ple->error_message,
                       ple->error_message); /* DLM */
            
            class_call(lensing_d3m2(mu,num_mu,ple->l_unlensed_max,d3m2),
                       ple->error_message,
                       ple->error_message); /* DLM */
            
        }
        
        class_call(lensing_d40(mu,num_mu,ple->l_unlensed_max,d40),
                   ple->error_message,
                   ple->error_message);
        
        class_call(lensing_d4m4(mu,num_mu,ple->l_unlensed_max,d4m4),
                   ple->error_message,
                   ple->error_message);
    }
    
    if (ple->calculate_pderivaties == _TRUE_ && ple->has_delensed_cls == _TRUE_){
        
        /** - Allocate secondary contiguous buffer **/
        class_alloc(buf2_dxx,
                    (ple->dl_size+1)*(ple->dl_size+1)*4*sizeof(double), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(buf4_dxx,
                    (ple->dl_size+1)*(ple->dl_size+1)*4*sizeof(double), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(buf3_dxx,
                    (ple->dl_size+1)*(ple->dl_size+1)*4*sizeof(double), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ple->cl_dl_tt_pderv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ple->cl_dl_te_pderv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ple->cl_dl_ee_pderv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ple->cl_dl_bb_pderv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ksi_dl_derv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ksiX_dl_derv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ksip_dl_derv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ksim_dl_derv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        
        class_alloc(ple->ddcl_dl_tt_pderv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ple->ddcl_dl_te_pderv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ple->ddcl_dl_ee_pderv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        class_alloc(ple->ddcl_dl_bb_pderv,
                    ple->dl_size*sizeof(double*), //make this precise! 5->4 and adjustments below.
                    ple->error_message);
        
        for(index_ll=0;index_ll<ple->dl_size;index_ll++){
            ple->cl_dl_tt_pderv[index_ll] = &(buf2_dxx[index_ll                   * (ple->dl_size+1)]);
            ple->cl_dl_te_pderv[index_ll] = &(buf2_dxx[(index_ll+ple->dl_size)   * (ple->dl_size+1)]);
            ple->cl_dl_ee_pderv[index_ll] = &(buf2_dxx[(index_ll+2*ple->dl_size) * (ple->dl_size+1)]);
            ple->cl_dl_bb_pderv[index_ll] = &(buf2_dxx[(index_ll+3*ple->dl_size) * (ple->dl_size+1)]);
        }
        for(index_ll=0;index_ll<ple->dl_size;index_ll++){
            ple->ddcl_dl_tt_pderv[index_ll] = &(buf4_dxx[index_ll                     * (ple->dl_size+1)]);
            ple->ddcl_dl_te_pderv[index_ll] = &(buf4_dxx[(index_ll+ple->dl_size)   * (ple->dl_size+1)]);
            ple->ddcl_dl_ee_pderv[index_ll] = &(buf4_dxx[(index_ll+2*ple->dl_size) * (ple->dl_size+1)]);
            ple->ddcl_dl_bb_pderv[index_ll] = &(buf4_dxx[(index_ll+3*ple->dl_size) * (ple->dl_size+1)]);
        }
        
        for(index_ll=0;index_ll<ple->dl_size;index_ll++){
            
            class_calloc(ksi_dl_derv[index_ll],
                         (num_mu-1),
                         sizeof(double),
                         ple->error_message);
            
            class_calloc(ksiX_dl_derv[index_ll],
                         (num_mu-1),
                         sizeof(double),
                         ple->error_message);
            
            class_calloc(ksip_dl_derv[index_ll],
                         (num_mu-1),
                         sizeof(double),
                         ple->error_message);
            
            class_calloc(ksim_dl_derv[index_ll],
                         (num_mu-1),
                         sizeof(double),
                         ple->error_message);
        }
    }
    
    if(ple->calculate_derviaties_wrt_unlensed == _TRUE_){
        if(ple->lensed_wrt_unlensed == _TRUE_){
            
            class_alloc(ksi_ln_derv,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ksiX_ln_derv,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ksip_ln_dervE,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ksip_ln_dervB,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ksim_ln_dervE,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ksim_ln_dervB,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            
            class_alloc(ple->cl_lens_derv_TT_TT,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_TE_TE,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_EE_EE,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_EE_BB,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_BB_EE,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_BB_BB,
                        (ple->l_unlensed_max+1)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            
            /*
            class_alloc(ple->cl_lens_derv_TT_TT,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_TE_TE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_EE_EE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_EE_BB,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_BB_EE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->cl_lens_derv_BB_BB,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            
            class_alloc(ple->ddcl_lens_derv_TT_TT,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_lens_derv_TE_TE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_lens_derv_EE_EE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_lens_derv_EE_BB,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_lens_derv_BB_EE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_lens_derv_BB_BB,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            */
        
        
        }        
        if(ple->delensed_wrt_unlensed == _TRUE_){
            if (ple->delensing_verbose > 2)  printf("the delensed_wrt_unlensed is true!\n"); /* DLM */
            
            class_alloc(ksi_dlu_derv,
                        (ple->l_unlensed_max+1)*sizeof(double*), 
                        ple->error_message);
            class_alloc(ksiX_dlu_derv,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            class_alloc(ksip_dlu_dervE,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            class_alloc(ksip_dlu_dervB,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            class_alloc(ksim_dlu_dervB,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            class_alloc(ksim_dlu_dervE,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            
            class_alloc(ple->cl_delens_derv_TT_TT,
                        (ple->l_unlensed_max+1)*sizeof(double*), 
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_TE_TE,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_EE_EE,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_EE_BB,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_BB_EE,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_BB_BB,
                        (ple->l_unlensed_max+1)*sizeof(double*),
                        ple->error_message);
            
            /*
            class_alloc(ple->cl_delens_derv_TT_TT,
                        (ple->l_size)*sizeof(double*), 
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_TE_TE,
                        (ple->l_size)*sizeof(double*), 
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_EE_EE,
                        (ple->l_size)*sizeof(double*), 
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_EE_BB,
                        (ple->l_size)*sizeof(double*), 
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_BB_EE,
                        (ple->l_size)*sizeof(double*), 
                        ple->error_message);
            class_alloc(ple->cl_delens_derv_BB_BB,
                        (ple->l_size)*sizeof(double*), 
                        ple->error_message);
            
            class_alloc(ple->ddcl_delens_derv_TT_TT,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_delens_derv_TE_TE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_delens_derv_EE_EE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_delens_derv_EE_BB,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_delens_derv_BB_EE,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            class_alloc(ple->ddcl_delens_derv_BB_BB,
                        (ple->l_size)*sizeof(double*), //make this precise! 5->4 and adjustments below.
                        ple->error_message);
            */
        }
        
        if(ple->lensed_wrt_unlensed == _TRUE_){
            
            for(index_ll=0;index_ll<ple->l_unlensed_max+1;index_ll++){
                
                class_calloc(ksi_ln_derv[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksiX_ln_derv[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksip_ln_dervE[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksim_ln_dervE[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksip_ln_dervB[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksim_ln_dervB[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                // ====================
                class_calloc(ple->cl_lens_derv_TT_TT[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_TE_TE[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_EE_EE[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_EE_BB[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_BB_EE[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_BB_BB[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                // ================================
                
            }
        }
        if(ple->delensed_wrt_unlensed == _TRUE_){
            for(index_ll=0;index_ll<ple->l_unlensed_max+1;index_ll++){
                
                class_calloc(ksi_dlu_derv[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksiX_dlu_derv[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksip_dlu_dervE[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksim_dlu_dervE[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksim_dlu_dervB[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);

                class_calloc(ksip_dlu_dervB[index_ll],
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                             
                // ============================
                class_calloc(ple->cl_delens_derv_TT_TT[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_TE_TE[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_EE_EE[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_EE_BB[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_BB_EE[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_BB_BB[index_ll],
                             (ple->l_unlensed_max+1),
                             sizeof(double),
                             ple->error_message);
                 // ================================
                
            }
        }
        
        // ===================================
        /*
        if(ple->lensed_wrt_unlensed == _TRUE_){
            for(index_l=0; index_l<ple->l_size; index_l++){
                
                
                class_calloc(ple->cl_lens_derv_TT_TT[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_TE_TE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_EE_EE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_EE_BB[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_BB_EE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_lens_derv_BB_BB[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_lens_derv_TT_TT[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_lens_derv_TE_TE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_lens_derv_EE_EE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_lens_derv_EE_BB[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_lens_derv_BB_EE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_lens_derv_BB_BB[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
            }
        }
        
        // ===================================
        
        if(ple->delensed_wrt_unlensed == _TRUE_){
            for(index_l=0; index_l<ple->l_size; index_l++){
                
                class_calloc(ple->cl_delens_derv_TT_TT[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_TE_TE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_EE_EE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_EE_BB[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_BB_EE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->cl_delens_derv_BB_BB[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_delens_derv_TT_TT[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_delens_derv_TE_TE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_delens_derv_EE_EE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_delens_derv_EE_BB[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_delens_derv_BB_EE[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ple->ddcl_delens_derv_BB_BB[index_l],
                             (ple->l_size),
                             sizeof(double),
                             ple->error_message);
            }
        }
        
        */
        // ==============================
    }
    
    /** - compute \f$ Cgl(\mu)\f$, \f$ Cgl2(\mu) \f$ and sigma2(\f$\mu\f$) */
    class_alloc(Cgl,
                num_mu*sizeof(double),
                ple->error_message);
    
    class_alloc(Cgl2,
                num_mu*sizeof(double),
                ple->error_message);
    
    class_alloc(sigma2,
                (num_mu-1)*sizeof(double), /* Zero separation is omitted */
                ple->error_message);
    
    if (ple->has_delensed_cls == _TRUE_) {
        /** - DLM: compute \f$ Cgl2_obs(\mu)\f$, \f$ Cgl2_cross(\mu) \f$
         and sigma2_obs(\f$\mu\f$)  sigma2_cross(\f$\mu\f$) */
        
        class_alloc(Cgl_obs,
                    num_mu*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(Cgl_cross,
                    num_mu*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(Cgl2_obs,
                    num_mu*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(Cgl2_cross,
                    num_mu*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(sigma2_hhbar,
                    (num_mu-1)*sizeof(double), /* Zero separation is omitted */
                    ple->error_message); /* DLM */
        
        class_alloc(sigma2_hh,
                    (num_mu-1)*sizeof(double), /* Zero separation is omitted */
                    ple->error_message); /* DLM */
        
    }
    
    class_alloc(cl_unlensed,
                phr->ct_size*sizeof(double),
                ple->error_message);
    
    /** - Locally store unlensed temperature \f$ cl_{tt}\f$ and potential \f$ cl_{pp}\f$ spectra **/
    class_alloc(cl_tt,
                (ple->l_unlensed_max+1)*sizeof(double),
                ple->error_message);
    if (ple->has_te==_TRUE_) {
        class_alloc(cl_te,
                    (ple->l_unlensed_max+1)*sizeof(double),
                    ple->error_message);
    }
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
      cl_ee[l] = cl_unlensed[ple->index_lt_ee];
      cl_bb[l] = cl_unlensed[ple->index_lt_bb];
    }
  }

  for (index_md = 0; index_md < phr->md_size; index_md++) {

    if (phr->md_size > 1)
      free(cl_md[index_md]);

    if (phr->ic_size[index_md] > 1)
      free(cl_md_ic[index_md]);

  }

  free(cl_md_ic);
  free(cl_md);

  /** - Compute sigma2\f$(\mu)\f$ and Cgl2(\f$\mu\f$) **/

  //debut = omp_get_wtime();
#pragma omp parallel for                        \
  private (index_mu,l)                          \
  schedule (static)
  for (index_mu=0; index_mu<num_mu; index_mu++) {

    Cgl[index_mu]=0;
    Cgl2[index_mu]=0;

    for (l=2; l<=ple->l_unlensed_max; l++) {

      Cgl[index_mu] += (2.*l+1.)*l*(l+1.)*
        cl_pp[l]*d11[index_mu][l];

      Cgl2[index_mu] += (2.*l+1.)*l*(l+1.)*
        cl_pp[l]*d1m1[index_mu][l];

    }

    Cgl[index_mu] /= 4.*_PI_;
    Cgl2[index_mu] /= 4.*_PI_;

  }

  for (index_mu=0; index_mu<num_mu-1; index_mu++) {
    /* Cgl(1.0) - Cgl(mu) */
    sigma2[index_mu] = Cgl[num_mu-1] - Cgl[index_mu];
  }
  //fin = omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in Cgl,Cgl2,sigma2=%4.3f s\n",cpu_time);


  /** - compute ksi, ksi+, ksi-, ksiX */

  /** - --> ksi is for TT **/
  if (ple->has_tt==_TRUE_) {

    class_calloc(ksi,
                 (num_mu-1),
                 sizeof(double),
                 ple->error_message);
  }

  /** - --> ksiX is for TE **/
  if (ple->has_te==_TRUE_) {

    class_calloc(ksiX,
                 (num_mu-1),
                 sizeof(double),
                 ple->error_message);
  }

  /** - --> ksip, ksim for EE, BB **/
  if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

    class_calloc(ksip,
                 (num_mu-1),
                 sizeof(double),
                 ple->error_message);

    class_calloc(ksim,
                 (num_mu-1),
                 sizeof(double),
                 ple->error_message);
  }

  for (l=2;l<=ple->l_unlensed_max;l++) {

    ll = (double)l;
    sqrt1[l]=sqrt((ll+2)*(ll+1)*ll*(ll-1));
    sqrt2[l]=sqrt((ll+2)*(ll-1));
    sqrt3[l]=sqrt((ll+3)*(ll-2));
    sqrt4[l]=sqrt((ll+4)*(ll+3)*(ll-2.)*(ll-3));
    sqrt5[l]=sqrt(ll*(ll+1));
  }


  //debut = omp_get_wtime();
#pragma omp parallel for                                                \
  private (index_mu,l,ll,res,resX,resp,resm,lens,lensp,lensm,           \
           fac,fac1,X_000,X_p000,X_220,X_022,X_p022,X_121,X_132,X_242)	\
  schedule (static)

  for (index_mu=0;index_mu<num_mu-1;index_mu++) {

    for (l=2;l<=ple->l_unlensed_max;l++) {

      ll = (double)l;

      fac = ll*(ll+1)/4.;
      fac1 = (2*ll+1)/(4.*_PI_);

      /* In the following we will keep terms of the form (sigma2)^k*(Cgl2)^m
         with k+m <= 2 */

      X_000 = exp(-fac*sigma2[index_mu]);
      X_p000 = -fac*X_000;
      /* X_220 = 0.25*sqrt1[l] * exp(-(fac-0.5)*sigma2[index_mu]); */
      X_220 = 0.25*sqrt1[l] * X_000; /* Order 0 */
      /* next 5 lines useless, but avoid compiler warning 'may be used uninitialized' */
      X_242=0.;
      X_132=0.;
      X_121=0.;
      X_p022=0.;
      X_022=0.;

      if (ple->has_te==_TRUE_ || ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        /* X_022 = exp(-(fac-1.)*sigma2[index_mu]); */
        X_022 = X_000 * (1+sigma2[index_mu]*(1+0.5*sigma2[index_mu])); /* Order 2 */
        X_p022 = -(fac-1.)*X_022; /* Old versions were missing the
                                     minus sign in this line, which introduced a very small error
                                     on the high-l C_l^TE lensed spectrum [credits for bug fix:
                                     Selim Hotinli] */

        /* X_242 = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2[index_mu]); */
        X_242 = 0.25*sqrt4[l] * X_000; /* Order 0 */
        if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

          /* X_121 = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2[index_mu]);
             X_132 = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2[index_mu]); */
          X_121 = -0.5*sqrt2[l] * X_000 * (1+2./3.*sigma2[index_mu]); /* Order 1 */
          X_132 = -0.5*sqrt3[l] * X_000 * (1+5./3.*sigma2[index_mu]); /* Order 1 */
        }
      }


      if (ple->has_tt==_TRUE_) {

        res = fac1*cl_tt[l];

        lens = (X_000*X_000*d00[index_mu][l] +
                X_p000*X_p000*d1m1[index_mu][l]
                *Cgl2[index_mu]*8./(ll*(ll+1)) +
                (X_p000*X_p000*d00[index_mu][l] +
                 X_220*X_220*d2m2[index_mu][l])
                *Cgl2[index_mu]*Cgl2[index_mu]);
        if (ppr->accurate_lensing == _FALSE_) {
          /* Remove unlensed correlation function */
          lens -= d00[index_mu][l];
        }
        res *= lens;
        ksi[index_mu] += res;
      }

      if (ple->has_te==_TRUE_) {

        resX = fac1*cl_te[l];


        lens = ( X_022*X_000*d20[index_mu][l] +
                 Cgl2[index_mu]*2.*X_p000/sqrt5[l] *
                 (X_121*d11[index_mu][l] + X_132*d3m1[index_mu][l]) +
                 0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                 ( ( 2.*X_p022*X_p000+X_220*X_220 ) *
                   d20[index_mu][l] + X_220*X_242*d4m2[index_mu][l] ) );
        if (ppr->accurate_lensing == _FALSE_) {
          lens -= d20[index_mu][l];
        }
        resX *= lens;
        ksiX[index_mu] += resX;
      }

      if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {

        resp = fac1*(cl_ee[l]+cl_bb[l]);
        resm = fac1*(cl_ee[l]-cl_bb[l]);

        lensp = ( X_022*X_022*d22[index_mu][l] +
                  2.*Cgl2[index_mu]*X_132*X_121*d31[index_mu][l] +
                  Cgl2[index_mu]*Cgl2[index_mu] *
                  ( X_p022*X_p022*d22[index_mu][l] +
                    X_242*X_220*d40[index_mu][l] ) );

        lensm = ( X_022*X_022*d2m2[index_mu][l] +
                  Cgl2[index_mu] *
                  ( X_121*X_121*d1m1[index_mu][l] +
                    X_132*X_132*d3m3[index_mu][l] ) +
                  0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                  ( 2.*X_p022*X_p022*d2m2[index_mu][l] +
                    X_220*X_220*d00[index_mu][l] +
                    X_242*X_242*d4m4[index_mu][l] ) );
        if (ppr->accurate_lensing == _FALSE_) {
          lensp -= d22[index_mu][l];
          lensm -= d2m2[index_mu][l];
        }
        resp *= lensp;
        resm *= lensm;
        ksip[index_mu] += resp;
        ksim[index_mu] += resm;
      }
    }
  }
  //fin = omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in ksi=%4.3f s\n",cpu_time);


  /** - compute lensed \f$ C_l\f$'s by integration */
  //debut = omp_get_wtime();
  if (ple->has_tt==_TRUE_) {
    class_call(lensing_lensed_cl_tt(ksi,d00,w8,num_mu-1,ple),
               ple->error_message,
               ple->error_message);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_tt(ple,cl_tt),
                 ple->error_message,
                 ple->error_message);
    }
    
    class_alloc(cl_pp,
                (ple->l_unlensed_max+1)*sizeof(double),
                ple->error_message);
    
    class_alloc(cl_lensed,
                ple->lt_size*sizeof(double),
                ple->error_message);
    
    
    if (ple->has_delensed_cls == _TRUE_) { /* DLM */
        
        class_alloc(nl_lensed,
                    ple->nlt_size*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_delensed,
                    ple->dlt_size*sizeof(double),
                    ple->error_message); /* DLM */
        
        /*--adding the allocation function for the
         observed lensing potential cls--*/
        
        class_alloc(cl_pp_obs,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_pp_cross,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(gl,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(hl,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(hlbar,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_lens_tt,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(hl_P,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(hlbar_P,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_lens_te,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_lens_ee,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_lens_bb,
                    (ple->l_lensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_delens_tt,
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_delens_te,
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_delens_ee,
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cl_delens_bb,
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(cratios,
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        class_alloc(mvolds,
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        class_alloc(tempds,
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(min_varr_ln,                                 //* separate if's here */ << CHANGE THIS CONDITION
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(min_ebeb_ln,                                //* separate if's here */ << CHANGE THIS CONDITION
                    (ple->l_delensed_max+1)*sizeof(double),
                    ple->error_message); /* DLM */
        
    }
    
    /* DLM: Generating/reading the temperature noise spectra and lensing
     reconstuction noise spectra from the input. */
    if(ple->has_delensed_cls == _TRUE_){
        if(ple->temperature_noise_type == idealized_tn){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Calculating idealized cmb temperature noise spectra in the form: \n (N_l= square of the instrumental noise in mu K-radians ) * exp[ l * (l+1) * ( square of the beamsize in radians / (8 log 2))^2 ].\n"); /* DLM */
            
            class_call_except(idealized_temperature_noise_spectrum_init(ple,phr),
                              ple->error_message,
                              ple->error_message,
                              lensing_free(ple)); /* DLM */
            
        }
        else if(ple->temperature_noise_type == external_tn){ /* DLM */
            
            if (ple->delensing_verbose > 0)
                printf("Calculating cmb temperature noise spectra from external data file.\n");
            
            class_call_except(external_temperature_noise_spectrum_init(ple,phr),
                              ple->error_message,
                              ple->error_message,
                              lensing_free(ple)); /* DLM */
        }
        if(ple->polarization_noise_type == idealized_pn){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Calculating idealized cmb polarization noise spectra in the form: \n (N_l= square of the instrumental noise in mu K-radians ) * exp[ l * (l+1) * ( square of the beamsize in radians / (8 log 2))^2 ].\n"); /* DLM */
            
            class_call_except(idealized_polarization_noise_spectrum_init(ple,phr),
                              ple->error_message,
                              ple->error_message,
                              lensing_free(ple)); /* DLM */
        }
        else if(ple->polarization_noise_type == external_pn){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Calculating cmb polariztion noise spectra from external data file.\n");
            
            class_call_except(external_polarization_noise_spectrum_init(ple,phr),
                              ple->error_message,
                              ple->error_message,
                              lensing_free(ple)); /* DLM */
        }
        
        if(ple->lens_rec_noise_type == internal_rn){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Will calculate lensing reconstruction noise spectra estimate from lensed CMB spectra a la astro-ph/0301031.\n");
            
        }
        else if(ple->lens_rec_noise_type == external_rn){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Calculating lensing reconstruction noise spectra from external data file.\n");
            
            class_call_except(external_lens_recon_noise_spectrum_init(ple),
                              ple->error_message,
                              ple->error_message,
                              lensing_free(ple)); /* DLM */
        }
        if(ple->cmb_spectra_type == internal_cmb){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Using CMB spectra calculated by CLASS.\n");
            
        }
        else if(ple->cmb_spectra_type == external_cmb){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Using external CMB spectra from file.\n");
            
            class_call_except(external_cmb_spectra_init(ple,phr),
                              ple->error_message,
                              ple->error_message,
                              lensing_free(ple)); /* DLM */
        }
        if(ple->cmb_spectra_lensed_type == internal_cmb){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Using *lensed* CMB spectra calculated by CLASS.\n");
            
        }
        else if(ple->cmb_spectra_lensed_type == external_cmb){ /* DLM */
            
            if (ple->delensing_verbose > 0) /* DLM */
                printf("Using external *lensed* CMB spectra from file.\n");
            
            class_call_except(external_lensed_cmb_spectra_init(ple,phr),
                              ple->error_message,
                              ple->error_message,
                              lensing_free(ple)); /* DLM */
        }
    }
    
    class_alloc(cl_md_ic,
                phr->md_size*sizeof(double *),
                ple->error_message);
    
    class_alloc(cl_md,
                phr->md_size*sizeof(double *),
                ple->error_message);
    
    for (index_md = 0; index_md < phr->md_size; index_md++) {
        
        if (phr->md_size > 1)
            
            class_alloc(cl_md[index_md],
                        phr->ct_size*sizeof(double),
                        ple->error_message);
        
        if (phr->ic_size[index_md] > 1)
            
            class_alloc(cl_md_ic[index_md],
                        phr->ic_ic_size[index_md]*phr->ct_size*sizeof(double),
                        ple->error_message);
    }
    if(ple->cmb_spectra_type != external_cmb){
        for (l=2; l<=ple->l_unlensed_max; l++) {
            class_call(harmonic_cl_at_l(phr,l,cl_unlensed,cl_md,cl_md_ic),
                       phr->error_message,
                       ple->error_message);
            cl_tt[l] = cl_unlensed[ple->index_lt_tt];
            cl_pp[l] = cl_unlensed[ple->index_lt_pp];
            
            if (ple->has_te==_TRUE_) {
                cl_te[l] = cl_unlensed[ple->index_lt_te];
            }
            if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
                cl_ee[l] = cl_unlensed[ple->index_lt_ee];
                cl_bb[l] = cl_unlensed[ple->index_lt_bb];
            }
        }
    }
    else{ // get the unlensed from external data
        for (l=2; l<=ple->l_unlensed_max; l++) {
            cl_tt[l] = ple->cl_tt_ext[l];
            cl_pp[l] = ple->cl_pp_ext[l];
            
            /* printf("%i %lf \n", l, cl_tt[l]); */
            
            if (ple->has_te==_TRUE_) {
                cl_te[l] = ple->cl_te_ext[l];
            }
            if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
                cl_ee[l] = ple->cl_ee_ext[l];
                cl_bb[l] = ple->cl_bb_ext[l];
            }
        }
    }
    
    for (index_md = 0; index_md < phr->md_size; index_md++) {
        
        if (phr->md_size > 1)
            free(cl_md[index_md]);
        
        if (phr->ic_size[index_md] > 1)
            free(cl_md_ic[index_md]);
        
    }
    
    free(cl_md_ic);
    free(cl_md);
    
    
    if(ple->cmb_spectra_lensed_type != external_cmb){
        
        /** - Compute sigma2\f$(\mu)\f$ and Cgl2(\f$\mu\f$) **/
        
        //debut = omp_get_wtime();
#pragma omp parallel for                        \
private (index_mu,l)                          \
schedule (static)
        for (index_mu=0; index_mu<num_mu; index_mu++) {
            
            Cgl[index_mu]=0;
            Cgl2[index_mu]=0;
            
            for (l=2; l<=ple->l_unlensed_max; l++) {
                
                Cgl[index_mu] += (2.*l+1.)*l*(l+1.)*
                cl_pp[l]*d11[index_mu][l];
                
                Cgl2[index_mu] += (2.*l+1.)*l*(l+1.)*
                cl_pp[l]*d1m1[index_mu][l];
                
            }
            
            Cgl[index_mu] /= 4.*_PI_;
            Cgl2[index_mu] /= 4.*_PI_;
            
        }
        
        for (index_mu=0; index_mu<num_mu-1; index_mu++) {
            /* Cgl(1.0) - Cgl(mu) */
            sigma2[index_mu] = Cgl[num_mu-1] - Cgl[index_mu];
            
        }
        
        //fin = omp_get_wtime();
        //cpu_time = (fin-debut);
        //printf("time in Cgl,Cgl2,sigma2=%4.3f s\n",cpu_time);
        
        /** - compute ksi, ksi+, ksi-, ksiX */
        
        /** - --> ksi is for TT **/
        if (ple->has_tt==_TRUE_) {
            
            class_calloc(ksi,
                         (num_mu-1),
                         sizeof(double),
                         ple->error_message);
        }
        
        /** - --> ksiX is for TE **/
        if (ple->has_te==_TRUE_) {
            
            class_calloc(ksiX,
                         (num_mu-1),
                         sizeof(double),
                         ple->error_message);
        }
        
        /** - --> ksip, ksim for EE, BB **/
        if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
            
            class_calloc(ksip,
                         (num_mu-1),
                         sizeof(double),
                         ple->error_message);
            
            class_calloc(ksim,
                         (num_mu-1),
                         sizeof(double),
                         ple->error_message);
        }
        
        for (l=2;l<=ple->l_unlensed_max;l++) {
            
            ll = (double)l;
            sqrt1[l]=sqrt((ll+2)*(ll+1)*ll*(ll-1));
            sqrt2[l]=sqrt((ll+2)*(ll-1));
            sqrt3[l]=sqrt((ll+3)*(ll-2));
            sqrt4[l]=sqrt((ll+4)*(ll+3)*(ll-2.)*(ll-3));
            sqrt5[l]=sqrt(ll*(ll+1));
        }
        //debut = omp_get_wtime();
#pragma omp parallel for                                                \
private (index_mu,l,ll,res,resX,resp,resm,lens,lensp,lensm,           \
fac,fac1,X_000,X_p000,X_220,X_022,X_p022,X_121,X_132,X_242)    \
schedule (static)
        //         printf("Calculating the lensing stuff.\n");
        for (index_mu=0;index_mu<num_mu-1;index_mu++) {
            
            for (l=2;l<=ple->l_unlensed_max;l++) {
                
                ll = (double)l;
                
                fac = ll*(ll+1)/4.;
                fac1 = (2*ll+1)/(4.*_PI_);
                
                /* take parameters that depend only on l out for efficiency */
                
                /* In the following we will keep terms of the form (sigma2)^k*(Cgl2)^m
                 with k+m <= 2 */
                
                X_000 = exp(-fac*sigma2[index_mu]);
                X_p000 = -fac*X_000;
                
                X_220 = 0.25*sqrt1[l] * exp(-(fac-0.5)*sigma2[index_mu]);
                //      X_220 = 0.25*sqrt1[l] * X_000; /* Order 0 */
                /* next 5 lines useless, but avoid compiler warning 'may be used uninitialized' */
                X_242=0.;
                X_132=0.;
                X_121=0.;
                X_p022=0.;
                X_022=0.;
                
                if (ple->has_te==_TRUE_ || ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
                    X_022 = exp(-(fac-1.)*sigma2[index_mu]);
                    // X_022 = X_000 * (1+sigma2[index_mu]*(1+0.5*sigma2[index_mu])); /* Order 2 */
                    X_p022 = -(fac-1.)*X_022; /* Old versions were missing the
                                               minus sign in this line, which introduced a very small error
                                               on the high-l C_l^TE lensed spectrum [credits for bug fix:
                                               Selim Hotinli] */
                    
                    X_242 = 0.25*sqrt4[l]  * exp(-(fac-5./2.)*sigma2[index_mu]);
                    // X_242 = 0.25*sqrt4[l] * X_000; /* Order 0 */
                    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
                        
                        X_121 = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2[index_mu]);
                        X_132 = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2[index_mu]);
                        // X_121 = -0.5*sqrt2[l] * X_000 * (1+2./3.*sigma2[index_mu]); /* Order 1 */
                        // X_132 = -0.5*sqrt3[l] * X_000 * (1+5./3.*sigma2[index_mu]); /* Order 1 */
                    }
                }
                
                
                if (ple->has_tt==_TRUE_) {
                    
                    res = fac1*cl_tt[l];
                    
                    lens = (X_000*X_000*d00[index_mu][l] +
                            X_p000*X_p000*d1m1[index_mu][l]
                            *Cgl2[index_mu]*8./(ll*(ll+1)) +
                            (X_p000*X_p000*d00[index_mu][l] +
                             X_220*X_220*d2m2[index_mu][l])
                            *Cgl2[index_mu]*Cgl2[index_mu]);
                    if (ppr->accurate_lensing == _FALSE_) {
                        /* Remove unlensed correlation function */
                        lens -= d00[index_mu][l];
                    }
                    res *= lens;
                    ksi[index_mu] += res;
                    if (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->lensed_wrt_unlensed == _TRUE_){
                        // ksi_ln_derv[l][index_mu] = res/cl_tt[l];
                        ksi_ln_derv[l][index_mu] = fac1*lens;
                    }
                    
                }
                
                
                if (ple->has_te==_TRUE_) {
                    
                    resX = fac1*cl_te[l];
                    
                    lens = ( X_022*X_000*d20[index_mu][l] +
                            Cgl2[index_mu]*2.*X_p000/sqrt5[l] *
                            (X_121*d11[index_mu][l] + X_132*d3m1[index_mu][l]) +
                            0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                            ( ( 2.*X_p022*X_p000+X_220*X_220 ) *
                             d20[index_mu][l] + X_220*X_242*d4m2[index_mu][l] ) );
                    if (ppr->accurate_lensing == _FALSE_) {
                        lens -= d20[index_mu][l];
                    }
                    resX *= lens;
                    ksiX[index_mu] += resX;
                    if (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->lensed_wrt_unlensed == _TRUE_){
                        // ksiX_ln_derv[l][index_mu] = resX/cl_te[l];
                        ksiX_ln_derv[l][index_mu] = fac1*lens;
                    }
                }
                
                if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
                    
                    resp = fac1*(cl_ee[l]+cl_bb[l]);
                    resm = fac1*(cl_ee[l]-cl_bb[l]);
                    
                    lensp = ( X_022*X_022*d22[index_mu][l] +
                             2.*Cgl2[index_mu]*X_132*X_121*d31[index_mu][l] +
                             Cgl2[index_mu]*Cgl2[index_mu] *
                             ( X_p022*X_p022*d22[index_mu][l] +
                              X_242*X_220*d40[index_mu][l] ) );
                    
                    lensm = ( X_022*X_022*d2m2[index_mu][l] +
                             Cgl2[index_mu] *
                             ( X_121*X_121*d1m1[index_mu][l] +
                              X_132*X_132*d3m3[index_mu][l] ) +
                             0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                             ( 2.*X_p022*X_p022*d2m2[index_mu][l] +
                              X_220*X_220*d00[index_mu][l] +
                              X_242*X_242*d4m4[index_mu][l] ));
                    if (ppr->accurate_lensing == _FALSE_) {
                        lensp -= d22[index_mu][l];
                        lensm -= d2m2[index_mu][l];
                    }
                    resp *= lensp;
                    resm *= lensm;
                    ksip[index_mu] += resp;
                    ksim[index_mu] += resm;
                    
                    if(ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->lensed_wrt_unlensed == _TRUE_){
                        // ksip_ln_dervE[l][index_mu] = resp/(cl_ee[l]+cl_bb[l]);
                        // ksim_ln_dervE[l][index_mu] = resm/(cl_ee[l]-cl_bb[l]);
                        // 
                        // ksip_ln_dervB[l][index_mu] = resp/(cl_ee[l]+cl_bb[l]);
                        // ksim_ln_dervB[l][index_mu] =-resm/(cl_ee[l]-cl_bb[l]);
                        
                        ksip_ln_dervE[l][index_mu] = fac1*lensp;
                        ksim_ln_dervE[l][index_mu] = fac1*lensm;
                        
                        ksip_ln_dervB[l][index_mu] = fac1*lensp;
                        ksim_ln_dervB[l][index_mu] =-fac1*lensm;
                    }
                }
            }
        }
        if (ple->delensing_verbose > 2) /* DLM */ printf("Done calculating lensing.\n");
        
        //fin = omp_get_wtime();
        //cpu_time = (fin-debut);
        ///printf("time in ksi=%4.3f s\n",cpu_time);
        
        
        /** - compute lensed \f$ C_l\f$'s by integration */
        //debut = omp_get_wtime();
        if (ple->has_tt==_TRUE_) {
            class_call(lensing_lensed_cl_tt(ksi,d00,w8,num_mu-1,ple),
                       ple->error_message,
                       ple->error_message);
            if (ppr->accurate_lensing == _FALSE_) {
                class_call(lensing_addback_cl_tt(ple,cl_tt),
                           ple->error_message,
                           ple->error_message);
            }

            if (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->lensed_wrt_unlensed == _TRUE_){
                class_call(lensing_lensed_cl_tt_derv_all(ksi_ln_derv,d00,w8,num_mu-1,ple),
                          ple->error_message,
                          ple->error_message);
            }
        }
        
        if (ple->has_te==_TRUE_) {
            class_call(lensing_lensed_cl_te(ksiX,d20,w8,num_mu-1,ple),
                       ple->error_message,
                       ple->error_message);
            if (ppr->accurate_lensing == _FALSE_) {
                class_call(lensing_addback_cl_te(ple,cl_te),
                           ple->error_message,
                           ple->error_message);
            }
            if (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->lensed_wrt_unlensed == _TRUE_){
//                 class_call(lensing_lensed_cl_te_derv(ksiX_ln_derv,d20,w8,num_mu-1,ple),
//                           ple->error_message,
//                           ple->error_message);
                class_call(lensing_lensed_cl_te_derv_all(ksiX_ln_derv,d20,w8,num_mu-1,ple),
                           ple->error_message,
                           ple->error_message);
            }
        }
        
        if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
            
            class_call(lensing_lensed_cl_ee_bb(ksip,ksim,d22,d2m2,w8,num_mu-1,ple),
                       ple->error_message,
                       ple->error_message);
            if (ppr->accurate_lensing == _FALSE_) {
                class_call(lensing_addback_cl_ee_bb(ple,cl_ee,cl_bb),
                           ple->error_message,
                           ple->error_message);
            }
            
            if(ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->lensed_wrt_unlensed == _TRUE_){
                
//                 class_call(lensing_lensed_cl_ee_bb_dervE(ksip_ln_dervE,ksim_ln_dervE,d22,d2m2,w8,num_mu-1,ple),
//                            ple->error_message,
//                            ple->error_message);
                
//                 class_call(lensing_lensed_cl_ee_bb_dervB(ksip_ln_dervB,ksim_ln_dervB,d22,d2m2,w8,num_mu-1,ple),
//                            ple->error_message,
//                            ple->error_message);
               
                class_call(lensing_lensed_cl_ee_bb_dervE_all(ksip_ln_dervE,ksim_ln_dervE,d22,d2m2,w8,num_mu-1,ple),
                           ple->error_message,
                           ple->error_message);
                
                class_call(lensing_lensed_cl_ee_bb_dervB_all(ksip_ln_dervB,ksim_ln_dervB,d22,d2m2,w8,num_mu-1,ple),
                           ple->error_message,
                           ple->error_message);
            }
        }
        
        class_call(array_spline_table_lines(ple->l,
                                            ple->l_size,
                                            ple->cl_lens,
                                            ple->lt_size,
                                            ple->ddcl_lens,
                                            _SPLINE_EST_DERIV_,
                                            ple->error_message),
                   ple->error_message,
                   ple->error_message);
        
        //fin=omp_get_wtime();
        //cpu_time = (fin-debut);
        //printf("time in final lensing computation=%4.3f s\n",cpu_time);
        
        if (ple->has_delensed_cls == _TRUE_){ /* DLM */
            
            for (l=2; l<=ple->l_lensed_max; l++) {
                
                ll = (double)l;
                
                class_call(lensing_cl_at_l(ple,
                                           ll,
                                           cl_lensed),
                           ple->error_message,
                           ple->error_message);
                
                cl_lens_tt[l] = cl_lensed[ple->index_lt_tt];
                cl_lens_te[l] = cl_lensed[ple->index_lt_te];
                cl_lens_ee[l] = cl_lensed[ple->index_lt_ee];
                cl_lens_bb[l] = cl_lensed[ple->index_lt_bb];
                
                
                /*--Setting below initially the delensed spectra equal to lensed spectra--*/
                /*------------------------------------------------------------------------*/
                if (ple->has_dl_tt == _TRUE_) cl_delens_tt[l] = cl_lens_tt[l]; /* DLM */
                if (ple->has_dl_te == _TRUE_) cl_delens_te[l] = cl_lens_te[l]; /* DLM */
                if (ple->has_dl_bb == _TRUE_ || ple->has_dl_ee == _TRUE_){        /* DLM */
                    cl_delens_ee[l] = cl_lens_ee[l]; /* DLM */
                    cl_delens_bb[l] = cl_lens_bb[l]; /* DLM */
                }
            }
        }
    }
    
    else{ // get the lensed from external data
        if (ple->has_delensed_cls == _TRUE_){
            for (l=2; l<=ple->l_lensed_max; l++) {
                cl_lens_tt[l] = ple->cl_lensed_tt_ext[l];
                cl_pp[l] = ple->cl_lensed_pp_ext[l]; // note perhaps this redundant definition is not ideal
                
                if (ple->has_dl_te==_TRUE_) {
                    cl_lens_te[l] = ple->cl_lensed_te_ext[l];
                }
                if (ple->has_dl_ee==_TRUE_ || ple->has_dl_bb==_TRUE_) {
                    cl_lens_ee[l] = ple->cl_lensed_ee_ext[l];
                    cl_lens_bb[l] = ple->cl_lensed_bb_ext[l];
                }
                
                /*--Setting below initially the delensed spectra equal to lensed spectra--*/
                /*------------------------------------------------------------------------*/
                if (ple->has_dl_tt == _TRUE_) cl_delens_tt[l] = cl_lens_tt[l]; /* DLM */
                if (ple->has_dl_te == _TRUE_) cl_delens_te[l] = cl_lens_te[l]; /* DLM */
                if (ple->has_dl_bb == _TRUE_ || ple->has_dl_ee == _TRUE_){        /* DLM */
                    cl_delens_ee[l] = cl_lens_ee[l]; /* DLM */
                    cl_delens_bb[l] = cl_lens_bb[l]; /* DLM */
                }
            }
        }
    }
    /*----------------------------------------------------------------*/
    /*------delensing the CMB modes and lensing reconstruction--------*/
    /*----------------------------------------------------------------*/
    if (ple->has_delensed_cls == _TRUE_){ /* DLM */
        
        if (ple->delensing_verbose > 0) printf("Calculating the delensed spectra.\n");
        
        int index_l; /* DLM */
        
        /*----------------------- begin delensing ------------------------*/
        /*----------------------------------------------------------------*/
        for( itr_index = 0; itr_index < ple->max_itr_steps; itr_index++){ /* DLM */
            if (ple->convergence_type == total && itr_index > 1
                && cratio <= ple->convergence_criterion_itr)
            {
                if(print_conv == _TRUE_ && ple->delensing_verbose>0)
                {
                    printf("Minimum variance N_l spectrum has converged.\n");
                    print_conv = _FALSE_;
                }
                continue;
            }
            if(ple->convergence_type == every &&
               type_2_flag == _TRUE_ && itr_index > 1)
            {
                if(print_conv == _TRUE_ && ple->delensing_verbose>0)
                {
                    printf("Minimum variance N_l spectrum has converged.\n");
                    print_conv = _FALSE_;
                }
                continue;
            }
            
            if (ple->has_lens_noise_rcn == _FALSE_){ /* DLM */
                
                if(ple->lens_rec_noise_type == internal_rn){ /* DLM */
                    
                    if(ple->has_nl_eb_itr == _FALSE_ &&
                       ple->has_nl_altr_itr == _FALSE_ &&
                       (ple->has_nl_all == _TRUE_ || ple->has_nl_diag == _TRUE_)) /* DLM */
                    {
                        
                        lensing_reconstr_nl_minvar(
                                                   d3m3, d33,
                                                   d3m2, d32,
                                                   d3m1, d31,
                                                   d30,
                                                   d2m2, d22,
                                                   d2m1, d21,
                                                   d20,
                                                   d1m1, d11,
                                                   d01,
                                                   d00,
                                                   w8,
                                                   num_mu-1,
                                                   cl_delens_tt, cl_delens_tt,
                                                   cl_delens_te, cl_delens_te,
                                                   cl_delens_ee, cl_delens_ee,
                                                   cl_bb, cl_delens_bb,
                                                   ple,
                                                   phr
                                                   ); /* DLM */
                    }
                    else if(ple->has_nl_eb_itr == _FALSE_ &&
                            ple->has_nl_altr_itr == _TRUE_ &&
                            (ple->has_nl_all == _TRUE_ || ple->has_nl_diag == _TRUE_)){ /* DLM */
                        
                        
                        lensing_reconstr_nl_minvar(
                                                   d3m3, d33,
                                                   d3m2, d32,
                                                   d3m1, d31,
                                                   d30,
                                                   d2m2, d22,
                                                   d2m1, d21,
                                                   d20,
                                                   d1m1, d11,
                                                   d01,
                                                   d00,
                                                   w8,
                                                   num_mu-1,
                                                   cl_lens_tt,   cl_lens_tt,
                                                   cl_lens_te,   cl_lens_te,
                                                   cl_lens_ee,   cl_lens_ee,
                                                   cl_bb,         cl_delens_bb,
                                                   ple,
                                                   phr
                                                   );  /* DLM */
                    }
                    else if(ple->has_nl_eb == _TRUE_){
                        
                        lensing_reconstr_nl_EB(
                                               d3m3, d33,
                                               d3m2, d32,
                                               d3m1, d31,
                                               d2m2, d22,
                                               d2m1, d21,
                                               d1m1, d11,
                                               d00,
                                               w8,
                                               num_mu-1,
                                               cl_lens_ee, cl_lens_ee,
                                               cl_bb,          cl_delens_bb,
                                               ple,
                                               phr
                                               );  /* DLM */
                        
                        
                        class_call(array_spline_table_lines(ple->l_dl,
                                                            ple->dl_size,
                                                            ple->nl_rcn,
                                                            ple->nlt_size,
                                                            ple->ddnl_rcn,
                                                            _SPLINE_EST_DERIV_,
                                                            ple->error_message),
                                   ple->error_message,
                                   ple->error_message); /* DLM */
                        
                    }
                    else if(ple->has_nl_eb_itr == _TRUE_ &&
                            (ple->has_nl_all == _TRUE_ || ple->has_nl_diag == _TRUE_)){ /* DLM */
                        
                        
/* internal loop for iterating the EB-varriance weigthed lensing noise and BB-spectra estimate*/
                        
                        if ((ple->convergence_type == total && itr_index > 1
                             && cratio <= ple->convergence_criterion_itr) ||
                            (ple->convergence_type == every && itr_index > 1
                             && type_2_flag == _TRUE_)) continue;
                        
// DLM: calculating the EB-weigthed varriance term in the lensing noise estimate with (de)lensed bb-spectra at the 1st (2nd+) iterations.
                        lensing_reconstr_nl_EB(
                                               d3m3, d33,
                                               d3m2, d32,
                                               d3m1, d31,
                                               d2m2, d22,
                                               d2m1, d21,
                                               d1m1, d11,
                                               d00,
                                               w8,
                                               num_mu-1,
                                               cl_lens_ee, cl_lens_ee,
                                               cl_bb,        cl_delens_bb,
                                               ple,
                                               phr
                                               );  /* DLM */
                        
                        
                        class_call(array_spline_table_lines(ple->l_dl,
                                                            ple->dl_size,
                                                            ple->nl_rcn,
                                                            ple->nlt_size,
                                                            ple->ddnl_rcn,
                                                            _SPLINE_EST_DERIV_,
                                                            ple->error_message),
                                   ple->error_message,
                                   ple->error_message); /* DLM */
                        
                        cratio = 0.; /* DLM */
                        
                        if(itr_index>1) type_2_flag = _TRUE_;
                        
                        // DLM: test the convergence of the EB-weigthed varriance estimate on the go.
                        for (l=2; l<=ple->l_lensed_max; l++) { /* DLM */
                            
                            ll = (double)l;
                            
                            if(itr_index>0) {
                                mvoldEB      = min_ebeb_ln[l]; /* DLM */
                                mvolds[l] = min_ebeb_ln[l];
                            }
                            
                            class_call(lensing_reconst_nl_at_l(ple,
                                                               ll,
                                                               nl_lensed),
                                       ple->error_message,
                                       ple->error_message); /* DLM */
                            
                            min_ebeb_ln[l] = nl_lensed[ple->index_nl_eb]; /* DLM */
                            
                            if(itr_index>0 && index_l>0)
                            {
                                tempdEB = (min_ebeb_ln[l]-mvoldEB)/min_ebeb_ln[l];
                                cratio += sqrt(tempdEB*tempdEB);
                                
                                tempds[l]  = (min_ebeb_ln[l]-mvolds[l])/min_ebeb_ln[l];
                                cratios[l] = sqrt(tempds[l]*tempds[l]);
                                
                                if (ple->convergence_type == every && itr_index > 1) {
                                    
                                    if (cratios[l] >= ple->convergence_criterion_itr){
                                        type_2_flag = _FALSE_;
                                    }
                                }
                            }
                        }
                        
                        if ((ple->convergence_type == total && itr_index > 0
                             && cratio <= ple->convergence_criterion_itr) ||
                            (ple->convergence_type == every && itr_index > 0
                             && type_2_flag == _TRUE_))
                        {
                            
// DLM: once the iterated EB-varriance weighted term is converged, add this and the rest of the weigths to calculate the minimum varriance noise.
                            lensing_reconstr_nl_minvar(
                                                       d3m3, d33,
                                                       d3m2, d32,
                                                       d3m1, d31,
                                                       d30,
                                                       d2m2, d22,
                                                       d2m1, d21,
                                                       d20,
                                                       d1m1, d11,
                                                       d01,
                                                       d00,
                                                       w8,
                                                       num_mu-1,
                                                       cl_lens_tt, cl_lens_tt,
                                                       cl_lens_te, cl_lens_te,
                                                       cl_lens_ee, cl_lens_ee,
                                                       cl_bb,          cl_lens_bb,
                                                       ple,
                                                       phr
                                                       );  /* DLM */
                            
                            ple->has_nl_eb_itr = _FALSE_;
                            
                            itr_index = ple->max_itr_steps + 1;
                            
                        }
                    }
                    
                    if(ple->has_nl_all == _TRUE_ && ple->has_nl_eb_itr == _FALSE_){
                        
                        double TTTT,TTTE,TTEE,TEEE,EEEE,TETE,TBTB,TBEB,EBEB,BBBB; /* DLM */
                        
                        for (index_l=0; index_l < ple->dl_size; index_l++)
                        {
                            
                            TTTT = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_tt]; /* DLM */
                            EEEE = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_ee]; /* DLM */
                            TETE = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_te]; /* DLM */
                            TBTB = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_tb]; /* DLM */
                            EBEB = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_eb]; /* DLM */
                            BBBB = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_bb]; /* DLM */
                            TTTE = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_ttte]; /* DLM */
                            TTEE = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_ttee]; /* DLM */
                            TEEE = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_teee]; /* DLM */
                            TBEB = ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_tbeb]; /* DLM */
                            
                            ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_minvar]=1./(1./BBBB +
                                                                                        (pow(TBEB,2)*(pow(TEEE,2) +
                                                                                                      pow(TTEE - TTTE,2) +
                                                                                                      TETE*(2.*TTEE - TTTT) -
                                                                                                      2.*TEEE*(TTEE + TTTE - TTTT) -
                                                                                                      EEEE*(TETE - 2.*TTTE + TTTT)) +
                                                                                         2.*TBEB*(-2.*TEEE*TTEE*TTTE +
                                                                                                  EEEE*pow(TTTE,2) +
                                                                                                  pow(TEEE,2)*TTTT +
                                                                                                  TETE*(pow(TTEE,2) - EEEE*TTTT)) -
                                                                                         TBTB*(-2.*TEEE*TTEE*TTTE +
                                                                                               EEEE*pow(TTTE,2) +
                                                                                               pow(TEEE,2)*TTTT +
                                                                                               TETE*(pow(TTEE,2) -
                                                                                                     EEEE*TTTT)) -
                                                                                         EBEB*(TETE*pow(TTEE,2) -
                                                                                               2.*TEEE*TTEE*TTTE +
                                                                                               EEEE*pow(TTTE,2) +
                                                                                               pow(TEEE,2)*TTTT -
                                                                                               EEEE*TETE*TTTT +
                                                                                               TBTB*(pow(TEEE,2) +
                                                                                                     pow(TTEE - TTTE,2) +
                                                                                                     TETE*(2.*TTEE - TTTT) -
                                                                                                     2.*TEEE*(TTEE + TTTE - TTTT) -
                                                                                                     EEEE*(TETE - 2.*TTTE + TTTT))))/
                                                                                        ((pow(TBEB,2) - EBEB*TBTB)*
                                                                                         (-2.*TEEE*TTEE*TTTE + EEEE*pow(TTTE,2) +
                                                                                          pow(TEEE,2)*TTTT + TETE*(pow(TTEE,2) - EEEE*TTTT)))); /* DLM */
                            
                        }
                    }
                    else if(ple->has_nl_diag == _TRUE_ && ple->has_nl_eb_itr == _FALSE_){ /* DLM */
                        
                        for (index_l=0; index_l < ple->dl_size; index_l++)
                        {
                            
                            ll = (double)ple->l_dl[index_l];
                            
                            ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_minvar]=1./(
                                                                                        1./ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_tt]
                                                                                        +1./ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_te]
                                                                                        +1./ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_ee]
                                                                                        +1./ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_bb]
                                                                                        +1./ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_eb]
                                                                                        +1./ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_tb]
                                                                                        ); /* DLM */
                            
                            
                        }
                    }
                    else if(ple->has_nl_eb_itr == _TRUE_ || ple->has_nl_eb == _TRUE_ ){/* DLM */
                        
                        for (index_l=0; index_l < ple->dl_size; index_l++)
                        {
                            ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_minvar]=ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_eb];
                            
                            ple->has_lens_noise_rcn = _TRUE_; /* DLM */
                            
                        }
                    }
                    
                    ple->has_lens_noise_rcn = _TRUE_; /* DLM */
                }
            
                class_call(array_spline_table_lines(ple->l_dl,
                                                    ple->dl_size,
                                                    ple->nl_rcn,
                                                    ple->nlt_size,
                                                    ple->ddnl_rcn,
                                                    _SPLINE_EST_DERIV_,
                                                    ple->error_message),
                           ple->error_message,
                           ple->error_message); /* DLM */
                
                
                if (ple->has_nl_eb_itr == _FALSE_ && itr_index <= ple->max_itr_steps) cratio = 0.; /* DLM */
                
                for (l=2; l<=ple->l_lensed_max; l++) { /* DLM */
                    
                    ll = (double)l;
                    
                    if (ple->has_nl_eb_itr == _FALSE_
                        && itr_index <= ple->max_itr_steps){
                        
                        if(itr_index>0){
                            mvold = min_varr_ln[l]; /* DLM */
                            mvolds[l] = min_varr_ln[l];
                            
                        }
                    }
                    
                    class_call(lensing_reconst_nl_at_l(ple,
                                                       ll,
                                                       nl_lensed),
                               ple->error_message,
                               ple->error_message); /* DLM */
                    
                    min_varr_ln[l] = nl_lensed[ple->index_nl_minvar]; /* DLM */
                    
                    if (ple->has_nl_eb_itr == _FALSE_
                        && itr_index <= ple->max_itr_steps){
                        if(itr_index>0)
                        {
                            
                            tempd = (min_varr_ln[l]-mvold)/min_varr_ln[l];
                            cratio += sqrt(tempd*tempd);
                            
                            tempds[l]  = (min_varr_ln[l]-mvolds[l])/min_varr_ln[l];
                            cratios[l] = sqrt(tempds[l]*tempds[l]);
                            
                        }
                    }
                }
            }
            else
            {
                for (l=0; l<ple->l_lensed_max; l++) { /* DLM */
                    /* NEEDS IMPROVEMENT */ // << also would cause error below when setting up min_varr_ln
                    min_varr_ln[l] = ple->pk_rcn_ext[l]; /* DLM */
                    cratios[l] = 0;
                }
                for (index_l=0; index_l < ple->dl_size; index_l++)
                {
                    ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_minvar]=min_varr_ln[(int)ple->l_dl[index_l]];
                }
                cratio = 0;
            }
            

            
            for (l=2; l<=ple->l_lensed_max; l++) { /* DLM */
                
                if (ple->has_dl_pp == _TRUE_) {
                    
                    cl_pp_obs[l]   = cl_pp[l] + min_varr_ln[l];
                    cl_pp_cross[l] = cl_pp[l];
                    
                    gl[l] = cl_pp[l] / ( cl_pp[l] + min_varr_ln[l] );
                    
                }
            }
            if(ple->delensing_verbose > 2) {
                printf("[DEBUG] (Total) convergence ratio: multipole(l)-sum of (N_l[i+1]-N_l[i])/N_l[i]=%e for iteration i=%d\n",cratio,itr_index);
                if(itr_index==0) printf("[DEBUG] N_l is minimum-variance lensing noise reconstruction spectrum\n");
            }
            
            if (ple->has_dl_tt == _TRUE_ || ple->has_dl_te == _TRUE_ ||
                ple->has_dl_ee == _TRUE_) { /* DLM */
                
                for (index_mu=0; index_mu<num_mu; index_mu++) {
                    
                    Cgl_obs[index_mu]=0;
                    Cgl_cross[index_mu]=0;
                    Cgl2_obs[index_mu]=0;
                    Cgl2_cross[index_mu]=0;
                    
                    for (l=2; l<=ple->l_lensed_max; l++) {
                        
                        ll = (double)l;
                        
                        Cgl_obs[index_mu] += (2.*ll+1.)*ll*(ll+1.)*
                        gl[l]*gl[l]*cl_pp_obs[l]*d11[index_mu][l];
                        
                        Cgl2_obs[index_mu] += (2.*ll+1.)*ll*(ll+1.)*
                        gl[l]*gl[l]*cl_pp_obs[l]*d1m1[index_mu][l];
                        
                        Cgl_cross[index_mu] += (2.*ll+1.)*ll*(ll+1.)*
                        gl[l]*cl_pp_cross[l]*d11[index_mu][l];
                        
                        Cgl2_cross[index_mu] += (2.*ll+1.)*ll*(ll+1.)*
                        gl[l]*cl_pp_cross[l]*d1m1[index_mu][l];
                        
                    }
                    
                    Cgl_obs[index_mu]   /= 4.*_PI_;
                    Cgl2_obs[index_mu]  /= 4.*_PI_;
                    Cgl_cross[index_mu] /= 4.*_PI_;
                    Cgl2_cross[index_mu]/= 4.*_PI_;
                    
                }
            }
            
            /* Having calculated lensed spectra, we now calculate the filters */
            for (l=2; l<=ple->l_lensed_max; l++) { /* DLM */
                
                ll = (double)l;
                
                hl[l] = cl_lens_tt[l]/(cl_lens_tt[l] + ple->pk_tn[l]);
                
                hlbar[l] = (sqrt(1.-hl[l]*hl[l]*
                                 (1.-exp(-ll*(ll+1.)*
                                         Cgl_obs[num_mu-1]/2.)))
                            -hl[l]*exp(-ll*(ll+1.)*Cgl_obs[num_mu-1]/4.));
                
                
                
                if (ple->has_dl_te == _TRUE_ || ple->has_dl_ee == _TRUE_) {
                    
                    hl_P[l] = cl_lens_ee[l]/(cl_lens_ee[l] + ple->pk_pn[l]);
                    
                    hlbar_P[l] =(sqrt(1.-hl_P[l]*hl_P[l]*
                                      (1.-exp(-ll*(ll+1.)*
                                              Cgl_obs[num_mu-1]/2.)))
                                 -hl_P[l]*exp(-ll*(ll+1.)*Cgl_obs[num_mu-1]/4.));
                    
                    /* Testing more aggressive filtering - JRM 06-2020 */
                    /*
                     hl_P[l] = hl[l];
                     hlbar_P[l] = hlbar[l];
                     */
                }
                
            }
            
            for (index_mu=0; index_mu<num_mu-1; index_mu++) {
                
                sigma2_hhbar[index_mu] = sigma2[index_mu]-((Cgl_cross[num_mu-1]-Cgl_cross[index_mu])-
                                                           0.5*Cgl_obs[num_mu-1]);
                
                sigma2_hh[index_mu] = sigma2[index_mu]-(2.*(Cgl_cross[num_mu-1]-Cgl_cross[index_mu])-
                                                        (Cgl_obs[num_mu-1] - Cgl_obs[index_mu]));
                
            }
            
            if (ple->has_dl_tt==_TRUE_){
                
                class_calloc(ksi_dl,
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
            }
            
            /** - --> ksiX_dl is for TE **/
            if (ple->has_dl_te==_TRUE_) {
                
                class_calloc(ksiX_dl,
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
            }
            
            /** - --> ksip_dl, ksim_dl for EE, BB **/
            if (ple->has_dl_ee==_TRUE_ || ple->has_dl_bb==_TRUE_) {
                
                class_calloc(ksip_dl,
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
                class_calloc(ksim_dl,
                             (num_mu-1),
                             sizeof(double),
                             ple->error_message);
                
            }
            
            if (ple->has_dl_tt == _TRUE_ || ple->has_dl_te == _TRUE_ ||
                ple->has_dl_ee == _TRUE_) {
                
                
#pragma omp parallel for                        \
private (index_mu,l,ll,fac,fac1,X_242_hhbar,X_132_hhbar,                \
X_121_hhbar,X_p022_hhbar,X_022_hhbar,X_242_hh,                    \
X_132_hh,X_121_hh,X_p022_hh,X_022_hh,X_000,X_000_hhbar,        \
X_000_hh,X_p000,X_p000_hhbar,X_p000_hh,X_220,X_220_hhbar,        \
X_220_hh,X_022,X_p022,X_242,X_121,X_132,                        \
res_dl_ml,lens_dl,resX_dl,resp_dl,resm_dl,lensp_dl,lensm_dl,    \
facd, sigma2_derv,sigma2_hh_derv,sigma2_hhbar_derv,X_000_derv, \
X_000_hhbar_derv,X_000_hh_derv, X_p000_derv,X_p000_hhbar_derv, \
X_p000_hh_derv,X_022_derv,X_022_hhbar_derv,X_p022_hh_derv,        \
X_p022_derv,X_p022_hhbar_derv,X_022_hh_derv,X_220_derv,         \
X_220_hhbar_derv,X_220_hh_derv,X_121_derv,X_121_hhbar_derv,    \
X_121_hh_derv,X_132_derv,X_132_hhbar_derv,X_132_hh_derv,        \
X_242_derv,X_242_hhbar_derv,X_242_hh_derv,kk,k,index1_l,        \
res_dl_derv,resX_dl_derv,resm_dl_derv,resp_dl_derv,             \
lens_dl_derv,lensX_dl_derv,lensm_dl_derv,lensp_dl_derv)        \
schedule (static)
                for (index_mu=0;index_mu<num_mu-1;index_mu++) {
                    for (l=2;l<=ple->l_lensed_max;l++) {
                        
                        ll = (double)l;
                        
                        fac = ll*(ll+1)/4.;
                        fac1 = (2*ll+1)/(4.*_PI_);
                        
                        X_242_hhbar = 0.;
                        X_132_hhbar = 0.;
                        X_121_hhbar = 0.;
                        X_p022_hhbar= 0.;
                        X_022_hhbar = 0.;
                        
                        X_242_hh = 0.;
                        X_132_hh = 0.;
                        X_121_hh = 0.;
                        X_p022_hh= 0.;
                        X_022_hh = 0.;
                        
                        X_000         = exp(-fac*sigma2[index_mu]);
                        X_000_hhbar = exp(-fac*sigma2_hhbar[index_mu]);
                        X_000_hh     = exp(-fac*sigma2_hh[index_mu]);
                        
                        X_p000          = -fac*X_000;
                        X_p000_hhbar = -fac*X_000_hhbar;
                        X_p000_hh     = -fac*X_000_hh;
                        
                        X_220         = 0.25*sqrt1[l] * exp(-(fac-0.5)*sigma2[index_mu]);
                        X_220_hhbar = 0.25*sqrt1[l] * exp(-(fac-0.5)*sigma2_hhbar[index_mu]);
                        X_220_hh     = 0.25*sqrt1[l] * exp(-(fac-0.5)*sigma2_hh[index_mu]);
                        
                        X_022         = exp(-(fac-1.)*sigma2[index_mu]);
                        X_022_hhbar = exp(-(fac-1.)*sigma2_hhbar[index_mu]);
                        X_022_hh     = exp(-(fac-1.)*sigma2_hh[index_mu]);
                        
                        /*fixed CLASS bug below*/
                        X_p022          = -(fac-1.)*X_022;
                        X_p022_hhbar = -(fac-1.)*X_022_hhbar;
                        X_p022_hh      = -(fac-1.)*X_022_hh;
                        
                        X_242         = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2[index_mu]);
                        X_242_hhbar = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2_hhbar[index_mu]);
                        X_242_hh    = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2_hh[index_mu]);
                        
                        if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
                            
                            X_121   = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2[index_mu]);
                            X_121_hhbar = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2_hhbar[index_mu]);
                            X_121_hh    = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2_hh[index_mu]);
                            
                            X_132         = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2[index_mu]);
                            X_132_hhbar = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2_hhbar[index_mu]);  // JM 2018 fixed typo here //
                            X_132_hh    = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2_hh[index_mu]);
                            
                        }
                        
                        if (ple->has_dl_tt == _TRUE_) {
                            
                            res_dl_ml = fac1*cl_tt[l];
                            
                            lens_dl = 0.;
                            
                            lens_dl += hlbar[l]*hlbar[l]
                            *( X_000*X_000*d00[index_mu][l]
                              +X_p000*X_p000*d1m1[index_mu][l]
                              *Cgl2[index_mu]*8./(ll*(ll+1))
                              +( X_p000*X_p000*d00[index_mu][l]
                                +X_220*X_220*d2m2[index_mu][l])
                              *Cgl2[index_mu]*Cgl2[index_mu]);
                            
                            lens_dl += 2.*hlbar[l]*hl[l]
                            *( X_000_hhbar*X_000_hhbar*d00[index_mu][l]
                              +(Cgl2[index_mu]-Cgl2_cross[index_mu])
                              *X_p000_hhbar*X_p000_hhbar*d1m1[index_mu][l]*8./(ll*(ll+1.))
                              +(Cgl2[index_mu]-Cgl2_cross[index_mu])
                              *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                              *(X_p000_hhbar*X_p000_hhbar*d00[index_mu][l]
                                +X_220_hhbar*X_220_hhbar*d2m2[index_mu][l]));
                            
                            lens_dl += hl[l]*hl[l]
                            *(X_000_hh*X_000_hh*d00[index_mu][l]+
                              (Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+Cgl2_obs[index_mu])
                              *X_p000_hh*X_p000_hh*d1m1[index_mu][l]*8./(ll*(ll+1.))
                              +(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+Cgl2_obs[index_mu])
                              *(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+Cgl2_obs[index_mu])
                              *(X_p000_hh*X_p000_hh*d00[index_mu][l]
                                +X_220_hh*X_220_hh*d2m2[index_mu][l]));
                            
                            
                            if (ppr->accurate_lensing == _FALSE_) {
                                /* Remove unlensed correlation function */
                                lens_dl -= d00[index_mu][l];
                            }
                            
                            res_dl_ml *= lens_dl;
                            
                            ksi_dl[index_mu] += res_dl_ml;
                            
                            if(ple->calculate_derviaties_wrt_unlensed == _TRUE_
                               && ple->delensed_wrt_unlensed == _TRUE_){
                                // ksi_dlu_derv[l][index_mu] = res_dl_ml/cl_tt[l];
                                ksi_dlu_derv[l][index_mu] = fac1*lens_dl;
                            }
                            
                        }
                        
                        if (ple->has_dl_te == _TRUE_) {
                            
                            resX_dl = fac1*cl_te[l];
                            
                            lens_dl = 0.;
                            
                            /*JRM: TE filtering; Changed from
                             lens_dl += hlbar_P[l]*hlbar[l]* */
                            lens_dl += hlbar_P[l]*hlbar_P[l]*
                            ( X_022*X_000*d20[index_mu][l] +
                             Cgl2[index_mu]*2.*X_p000/sqrt5[l] *
                             (X_121*d11[index_mu][l] + X_132*d3m1[index_mu][l]) +
                             0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                             ( ( 2.*X_p022*X_p000+X_220*X_220 ) *
                              d20[index_mu][l] + X_220*X_242*d4m2[index_mu][l]));
                            
                            
                            /*JRM: TE filtering; Changed from
                             lens_dl += (hlbar_P[l]*hl[l]+hlbar[l]*hl_P[l])* */
                            lens_dl += (hlbar_P[l]*hl_P[l]+hlbar_P[l]*hl_P[l])*
                            (X_022_hhbar*X_000_hhbar*d20[index_mu][l]+
                             (Cgl2[index_mu]-Cgl2_cross[index_mu])*2./sqrt5[l]*
                             X_p000_hhbar*(X_121_hhbar*d11[index_mu][l]+
                                           X_132_hhbar*d3m1[index_mu][l])+
                             0.5*(Cgl2[index_mu]-Cgl2_cross[index_mu])*
                             (Cgl2[index_mu]-Cgl2_cross[index_mu])*
                             ((2.*X_p022_hhbar*X_p000_hhbar+
                               X_220_hhbar*X_220_hhbar)*d20[index_mu][l]+
                              X_220_hhbar*X_242_hhbar*d4m2[index_mu][l]));
                            
                            /*JRM: TE filtering; Changed from
                             lens_dl += hl_P[l]*hl[l]* */
                            lens_dl += hl_P[l]*hl_P[l]*
                            (X_022_hh*X_000_hh*d20[index_mu][l]+
                             (Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+
                              Cgl2_obs[index_mu])*2./sqrt5[l]*
                             X_p000_hh*(X_121_hh*d11[index_mu][l]+
                                        X_132_hh*d3m1[index_mu][l])+
                             0.5*(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+
                                  Cgl2_obs[index_mu])*
                             (Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+
                              Cgl2_obs[index_mu])*
                             ((2.*X_p022_hh*X_p000_hh+
                               X_220_hh*X_220_hh)*d20[index_mu][l]+
                              X_220_hh*X_242_hh*d4m2[index_mu][l]));
                            
                            if (ppr->accurate_lensing == _FALSE_) {
                                lens_dl -= d20[index_mu][l];
                            }
                            
                            resX_dl *= lens_dl;
                            
                            ksiX_dl[index_mu] += resX_dl;
                            
                            if(ple->calculate_derviaties_wrt_unlensed == _TRUE_
                               && ple->delensed_wrt_unlensed == _TRUE_){
                                // ksiX_dlu_derv[l][index_mu] = resX_dl/cl_te[l];
                                ksiX_dlu_derv[l][index_mu] = fac1*lens_dl;
                            }
                            
                        }
                        
                        if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
                            
                            resp_dl = fac1*(cl_ee[l]+cl_bb[l]);
                            resm_dl = fac1*(cl_ee[l]-cl_bb[l]);
                            
                            lensp_dl = 0.;
                            
                            lensp_dl += hlbar_P[l]*hlbar_P[l]*
                            ( X_022*X_022*d22[index_mu][l] +
                             2.*Cgl2[index_mu]*X_132*X_121*d31[index_mu][l] +
                             Cgl2[index_mu]*Cgl2[index_mu] *
                             ( X_p022*X_p022*d22[index_mu][l] +
                              X_242*X_220*d40[index_mu][l] ) );
                            
                            lensp_dl += 2.*hlbar_P[l]*hl_P[l]*
                            (X_022_hhbar*X_022_hhbar*d22[index_mu][l]+
                             2.*(Cgl2[index_mu]-Cgl2_cross[index_mu])*
                             X_132_hhbar*X_121_hhbar*d31[index_mu][l]+
                             (Cgl2[index_mu]-Cgl2_cross[index_mu])*
                             (Cgl2[index_mu]-Cgl2_cross[index_mu])*
                             (X_p022_hhbar*X_p022_hhbar*d22[index_mu][l]+
                              X_242_hhbar*X_220_hhbar*d40[index_mu][l]));
                            
                            lensp_dl += hl_P[l]*hl_P[l]*
                            (X_022_hh*X_022_hh*d22[index_mu][l]+
                             2.*(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+
                                 Cgl2_obs[index_mu])*
                             X_132_hh*X_121_hh*d31[index_mu][l]+
                             (Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+
                              Cgl2_obs[index_mu])*
                             (Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+
                              Cgl2_obs[index_mu])*
                             (X_p022_hh*X_p022_hh*d22[index_mu][l]+
                              X_242_hh*X_220_hh*d40[index_mu][l]));
                            
                            
                            if (ppr->accurate_lensing == _FALSE_)
                            {
                                lensp_dl -= d22[index_mu][l];
                            }
                            
                            resp_dl *= lensp_dl;
                            
                            ksip_dl[index_mu] += resp_dl;
                            
                            lensm_dl = 0.;
                            
                            lensm_dl += hlbar_P[l]*hlbar_P[l]
                            *( X_022*X_022*d2m2[index_mu][l]
                              +Cgl2[index_mu]
                              *( X_121*X_121*d1m1[index_mu][l]
                                +X_132*X_132*d3m3[index_mu][l])
                              +0.5*Cgl2[index_mu]*Cgl2[index_mu]
                              *(2.*X_p022*X_p022*d2m2[index_mu][l]
                                +X_220*X_220*d00[index_mu][l]
                                +X_242*X_242*d4m4[index_mu][l]));
                            
                            lensm_dl += 2.*hlbar_P[l]*hl_P[l]
                            *(X_022_hhbar*X_022_hhbar*d2m2[index_mu][l]
                              +(Cgl2[index_mu]-Cgl2_cross[index_mu])
                              *(X_121_hhbar*X_121_hhbar*d1m1[index_mu][l]
                                +X_132_hhbar*X_132_hhbar*d3m3[index_mu][l])
                              +0.5*(Cgl2[index_mu]-Cgl2_cross[index_mu])
                              *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                              *(2.*X_p022_hhbar*X_p022_hhbar*d2m2[index_mu][l]
                                +X_220_hhbar*X_220_hhbar*d00[index_mu][l]
                                +X_242_hhbar*X_242_hhbar*d4m4[index_mu][l]));
                            
                            lensm_dl += hl_P[l]*hl_P[l]*
                            (X_022_hh*X_022_hh*d2m2[index_mu][l]
                             +(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                               +Cgl2_obs[index_mu])
                             *(X_121_hh*X_121_hh*d1m1[index_mu][l]
                               +X_132_hh*X_132_hh*d3m3[index_mu][l])
                             +0.5*( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                   +Cgl2_obs[index_mu])
                             *(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                               +Cgl2_obs[index_mu])
                             *(2.*X_p022_hh*X_p022_hh*d2m2[index_mu][l]
                               +X_220_hh*X_220_hh*d00[index_mu][l]
                               +X_242_hh*X_242_hh*d4m4[index_mu][l]));
                            
                            resm_dl *= lensm_dl;
                            
                            ksim_dl[index_mu] += resm_dl;
                            
                            if(ple->calculate_derviaties_wrt_unlensed == _TRUE_
                               && ple->delensed_wrt_unlensed == _TRUE_){
                                // ksim_dlu_dervE[l][index_mu] = resm_dl/(cl_ee[l]-cl_bb[l]);
                                // ksim_dlu_dervB[l][index_mu] = -resm_dl/(cl_ee[l]-cl_bb[l]);
                                
                                // ksip_dlu_dervE[l][index_mu] = resp_dl/(cl_ee[l]+cl_bb[l]);
                                // ksip_dlu_dervB[l][index_mu] = resp_dl/(cl_ee[l]+cl_bb[l]);
                                
                                ksim_dlu_dervE[l][index_mu] = fac1*lensm_dl;
                                ksim_dlu_dervB[l][index_mu] = -fac1*lensm_dl;
                                
                                ksip_dlu_dervE[l][index_mu] = fac1*lensp_dl;
                                ksip_dlu_dervB[l][index_mu] = fac1*lensp_dl;
                            }
                            
                        }
                        
                        if(((ple->convergence_type == total && itr_index > 1
                             && cratio <= ple->convergence_criterion_itr) ||
                            (ple->convergence_type == every && itr_index > 1
                             && type_2_flag == _TRUE_) ||
                            ple->lens_rec_noise_type == external_rn || 
                           ple->has_itr_delensing == _FALSE_) &&
                           ple->calculate_pderivaties == _TRUE_)
                        {
                            
                            for(index1_l=0; index1_l<ple->dl_size; index1_l++){
                                
                                k = (int)ple->l_dl[index1_l];
                                
                                kk = (double)ple->l_dl[index1_l];
                                
                                facd = (2.*kk+1.)*(kk+1.)*kk/(4.*_PI_);
                                
                                sigma2_derv       = facd*( d11[num_mu-1][k]
                                                          -d11[index_mu][k]);
                                
                                sigma2_hhbar_derv = sigma2_derv-facd*(gl[k]*( d11[num_mu-1][k]
                                                                             -d11[index_mu][k])
                                                                      -0.5*gl[k]*gl[k]*d11[num_mu-1][k]);
                                
                                sigma2_hh_derv       = sigma2_derv-facd*(2.*gl[k]*( d11[num_mu-1][k]
                                                                                   -d11[index_mu][k])
                                                                         -gl[k]*gl[k]*( d11[num_mu-1][k]
                                                                                       -d11[index_mu][k]));
                                
                                X_000_derv=0.;  X_p000_derv=0.; X_220_derv=0.; X_022_derv=0.;
                                X_p022_derv=0.; X_121_derv=0.;  X_132_derv=0.; X_242_derv=0.;
                                
                                X_000_hh_derv=0.;  X_p000_hh_derv=0.; X_220_hh_derv=0.; X_022_hh_derv=0.;
                                X_p022_hh_derv=0.; X_121_hh_derv=0.;  X_132_hh_derv=0.; X_242_hh_derv=0.;
                                
                                X_000_hhbar_derv=0.;  X_p000_hhbar_derv=0.; X_220_hhbar_derv=0.; X_022_hhbar_derv=0.;
                                X_p022_hhbar_derv=0.; X_121_hhbar_derv=0.;  X_132_hhbar_derv=0.; X_242_hhbar_derv=0.;
                                
                                X_000_derv           = -fac*sigma2_derv*X_000;
                                X_000_hhbar_derv  = -fac*sigma2_hhbar_derv*X_000_hhbar;
                                X_000_hh_derv      = -fac*sigma2_hh_derv*X_000_hh;
                                
                                X_p000_derv       = -fac*X_000_derv;
                                X_p000_hhbar_derv = -fac*X_000_hhbar_derv;
                                X_p000_hh_derv       = -fac*X_000_hh_derv;
                                
                                X_220_derv          = -(fac-0.5)*sigma2_derv*X_220;
                                X_220_hhbar_derv = -(fac-0.5)*sigma2_hhbar_derv*X_220_hhbar;
                                X_220_hh_derv      = -(fac-0.5)*sigma2_hh_derv*X_220_hh;
                                
                                X_022_derv             = -(fac-1.)*sigma2_derv*X_022;
                                X_022_hhbar_derv  = -(fac-1.)*sigma2_hhbar_derv*X_022_hhbar;
                                X_022_hh_derv       = -(fac-1.)*sigma2_hh_derv*X_022_hh;
                                
                                X_p022_derv       = -(fac-1.)*X_022_derv;
                                X_p022_hhbar_derv = -(fac-1.)*X_022_hhbar_derv;
                                X_p022_hh_derv       = -(fac-1.)*X_022_hh_derv;
                                
                                X_121_derv            = -(fac-2./3.)*sigma2_derv*X_121;
                                X_121_hhbar_derv = -(fac-2./3.)*sigma2_hhbar_derv*X_121_hhbar;
                                X_121_hh_derv      = -(fac-2./3.)*sigma2_hh_derv*X_121_hh;
                                
                                X_132_derv           = -(fac-5./3.)*sigma2_derv*X_132;
                                X_132_hhbar_derv = -(fac-5./3.)*sigma2_hhbar_derv*X_132_hhbar;
                                X_132_hh_derv      = -(fac-5./3.)*sigma2_hh_derv*X_132_hh;
                                
                                X_242_derv            = -(fac-5./2.)*sigma2_derv*X_242;
                                X_242_hhbar_derv = -(fac-5./2.)*sigma2_hhbar_derv*X_242_hhbar;
                                X_242_hh_derv      = -(fac-5./2.)*sigma2_hh_derv*X_242_hh;
                                
                                /*------------ Cl_TT ------------*/
                                
                                res_dl_derv = fac1*cl_tt[l];
                                
                                lens_dl_derv = 0.;
                                
                                lens_dl_derv += hlbar[l]*hlbar[l]
                                *(2.*X_000*X_000_derv*d00[index_mu][l]                        //1
                                  +2.*X_p000*X_p000_derv*d1m1[index_mu][l]                    //2
                                  *Cgl2[index_mu]*8./(ll*(ll+1.))                            //2
                                  +X_p000*X_p000*d1m1[index_mu][l]                            //2
                                  *facd*d1m1[index_mu][k]*8./(ll*(ll+1.))                    //2
                                  +( 2.*X_p000*X_p000_derv*d00[index_mu][l]                 //3
                                    +2.*X_220*X_220_derv*d2m2[index_mu][l])                 //3
                                  *Cgl2[index_mu]*Cgl2[index_mu]                            //3
                                  +( X_p000*X_p000*d00[index_mu][l]                            //3
                                    +X_220*X_220*d2m2[index_mu][l])                            //3
                                  *2.*facd*d1m1[index_mu][k]*Cgl2[index_mu]);
                                
                                if(ple->derv_type == lensed) lens_dl_derv /= hlbar[l]*hlbar[l];
                                
                                else if(ple->derv_type == delensed)
                                {
                                    
                                    lens_dl_derv += 2.*hlbar[l]*hl[l]
                                    *(2.*X_000_hhbar*X_000_hhbar_derv*d00[index_mu][l]
                                      +(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                      *2.*X_p000_hhbar*X_p000_hhbar_derv*d1m1[index_mu][l]*8./(ll*(ll+1.))
                                      +(facd*(1.-gl[k])*d1m1[index_mu][k])
                                      *X_p000_hhbar*X_p000_hhbar*d1m1[index_mu][l]*8./(ll*(ll+1.)) // silly bug here!
                                      +( 2.*X_p000_hhbar*X_p000_hhbar_derv*d00[index_mu][l]
                                        +2.*X_220_hhbar*X_220_hhbar_derv*d2m2[index_mu][l])
                                      *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                      *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                      +2.*(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                      *(facd*(1.-gl[k])*d1m1[index_mu][k]) // d(Cgl2-Cgl2_cross)/dCpp
                                      *( X_p000_hhbar*X_p000_hhbar*d00[index_mu][l]
                                        +X_220_hhbar*X_220_hhbar*d2m2[index_mu][l])
                                      );
                                    
                                    lens_dl_derv += hl[l]*hl[l]
                                    *( 2.*X_000_hh*X_000_hh_derv*d00[index_mu][l]
                                      +2.*X_p000_hh*X_p000_hh_derv*d1m1[index_mu][l]
                                      *(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                        +Cgl2_obs[index_mu])*8./(ll*(ll+1.))
                                      +X_p000_hh*X_p000_hh*d1m1[index_mu][l]*8./(ll*(ll+1.))
                                      *(facd*(1.-2.*gl[k]+gl[k]*gl[k])*d1m1[index_mu][k])
                                      +( 2.*X_p000_hh*X_p000_hh_derv*d00[index_mu][l]
                                        +2.*X_220_hh*X_220_hh_derv*d2m2[index_mu][l])
                                      *(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+Cgl2_obs[index_mu])
                                      *(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+Cgl2_obs[index_mu])
                                      +2.*(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]+Cgl2_obs[index_mu])
                                      *(facd*(1.-2.*gl[k]+gl[k]*gl[k])*d1m1[index_mu][k])
                                      *( X_p000_hh*X_p000_hh*d00[index_mu][l]
                                        +X_220_hh*X_220_hh*d2m2[index_mu][l])
                                      );
                                }
                                
                                /*------------ Cl_TE ------------*/
                                resX_dl_derv = fac1*cl_te[l];
                                
                                lensX_dl_derv = 0.;
                                
                                /*JRM: TE filtering; Changed from
                                 lensX_dl_derv += hlbar_P[l]*hlbar[l]* */
                                lensX_dl_derv += hlbar_P[l]*hlbar_P[l]*
                                (( X_022_derv*X_000
                                  +X_022*X_000_derv)*d20[index_mu][l]
                                 +2.*( facd*d1m1[index_mu][k]*X_p000
                                      +Cgl2[index_mu]*X_p000_derv)
                                 *( X_121*d11[index_mu][l]
                                   +X_132*d3m1[index_mu][l])/sqrt5[l]
                                 +2.*Cgl2[index_mu]*X_p000
                                 *( X_121_derv*d11[index_mu][l]
                                   +X_132_derv*d3m1[index_mu][l])/sqrt5[l]
                                 +(facd*d1m1[index_mu][k]*Cgl2[index_mu])
                                 *((2.*X_p022*X_p000+X_220*X_220)*d20[index_mu][l]
                                   +X_220*X_242*d4m2[index_mu][l])
                                 +0.5*Cgl2[index_mu]*Cgl2[index_mu]
                                 *(( 2.*X_p022_derv*X_p000
                                    +2.*X_p022*X_p000_derv
                                    +2.*X_220_derv*X_220)*d20[index_mu][l]
                                   +( X_220_derv*X_242
                                     +X_220*X_242_derv)*d4m2[index_mu][l]));
                                
                                
                                /*JRM: TE filtering; Changed from
                                 if(ple->derv_type == lensed) lensX_dl_derv /= hlbar_P[l]*hlbar[l]; */
                                if(ple->derv_type == lensed) lensX_dl_derv /= hlbar_P[l]*hlbar_P[l];
                                
                                else if(ple->derv_type == delensed)
                                {
                                    /*JRM: TE filtering; Changed from
                                     lensX_dl_derv += (hlbar_P[l]*hl[l]+hlbar[l]*hl_P[l])* */
                                    lensX_dl_derv += (hlbar_P[l]*hl_P[l]+hlbar_P[l]*hl_P[l])*
                                    (( X_022_hhbar_derv*X_000_hhbar
                                      +X_022_hhbar*X_000_hhbar_derv)*d20[index_mu][l]
                                     +(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *(( X_121_hhbar*d11[index_mu][l]
                                        +X_132_hhbar*d3m1[index_mu][l])
                                       *X_p000_hhbar_derv*2./sqrt5[l]
                                       +( X_121_hhbar_derv*d11[index_mu][l]
                                         +X_132_hhbar_derv*d3m1[index_mu][l])
                                       *X_p000_hhbar*2./sqrt5[l])
                                     +facd*(1.-gl[k])*d1m1[index_mu][k]
                                     *(( X_121_hhbar*d11[index_mu][l]
                                        +X_132_hhbar*d3m1[index_mu][l])
                                       *X_p000_hhbar*2./sqrt5[l])
                                     +(facd*(1.-gl[k])*d1m1[index_mu][k])
                                     *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *((2.*X_p022_hhbar*X_p000_hhbar+
                                        X_220_hhbar*X_220_hhbar)*d20[index_mu][l]
                                       +X_220_hhbar*X_242_hhbar*d4m2[index_mu][l])
                                     +0.5*(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *((2.*( X_p022_hhbar*X_p000_hhbar_derv
                                            +X_p022_hhbar_derv*X_p000_hhbar)
                                        +2.*X_220_hhbar*X_220_hhbar_derv)*d20[index_mu][l]
                                       +( X_220_hhbar*X_242_hhbar_derv
                                         +X_220_hhbar_derv*X_242_hhbar)*d4m2[index_mu][l]));
                                    
                                    
                                    /*JRM: TE filtering; Changed from
                                     lensX_dl_derv += hl_P[l]*hl[l]* */
                                    lensX_dl_derv += hl_P[l]*hl_P[l]*
                                    (( X_022_hh_derv*X_000_hh
                                      +X_022_hh*X_000_hh_derv)*d20[index_mu][l]
                                     +( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *(( X_121_hh*d11[index_mu][l]
                                        +X_132_hh*d3m1[index_mu][l])
                                       *X_p000_hh_derv
                                       +( X_121_hh_derv*d11[index_mu][l]
                                         +X_132_hh_derv*d3m1[index_mu][l])
                                       *X_p000_hh)*2./sqrt5[l]
                                     +(facd*(1.-2.*gl[k]+gl[k]*gl[k])
                                       *d1m1[index_mu][k])
                                     *( X_121_hh*d11[index_mu][l]
                                       +X_132_hh*d3m1[index_mu][l])
                                     *X_p000_hh*2./sqrt5[l]
                                     +0.5*( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                           +Cgl2_obs[index_mu])
                                     *( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *((2.*( X_p022_hh_derv*X_p000_hh
                                            +X_p022_hh*X_p000_hh_derv)
                                        +2.*X_220_hh_derv*X_220_hh)*d20[index_mu][l]
                                       +( X_220_hh_derv*X_242_hh
                                         +X_220_hh*X_242_hh_derv)*d4m2[index_mu][l])
                                     +(facd*(1.-2.*gl[k]+gl[k]*gl[k])*d1m1[index_mu][k])
                                     *( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *(( 2.*X_p022_hh*X_p000_hh
                                        +X_220_hh*X_220_hh)*d20[index_mu][l]
                                       +X_220_hh*X_242_hh*d4m2[index_mu][l]));
                                }
                                
                                
                                
                                /*------------ Cl_EE & CL_BB ------------*/
                                
                                resp_dl_derv = fac1*(cl_ee[l]+cl_bb[l]);
                                
                                lensp_dl_derv = 0.;
                                
                                lensp_dl_derv += hlbar_P[l]*hlbar_P[l]*
                                (2.*X_022*X_022_derv*d22[index_mu][l]
                                 +2.*( facd*d1m1[index_mu][k]
                                      *X_132*X_121*d31[index_mu][l])
                                 +2.*( Cgl2[index_mu]*d31[index_mu][l]
                                      *(X_132_derv*X_121+X_132*X_121_derv))
                                 +2.*facd*d1m1[index_mu][k]*Cgl2[index_mu]
                                 *( X_p022*X_p022*d22[index_mu][l]
                                   +X_242*X_220*d40[index_mu][l] )
                                 +2.*Cgl2[index_mu]*Cgl2[index_mu]
                                 *(2.*X_p022*X_p022_derv*d22[index_mu][l]
                                   +( X_242*X_220_derv
                                     +X_242_derv*X_220)*d40[index_mu][l]));
                                
                                
                                if(ple->derv_type == lensed) lensp_dl_derv /= hlbar_P[l]*hlbar_P[l];
                                
                                else if(ple->derv_type == delensed)
                                {
                                    
                                    lensp_dl_derv += 2.*hlbar_P[l]*hl_P[l]*
                                    ( 2.*X_022_hhbar*X_022_hhbar_derv*d22[index_mu][l]
                                     +2.*(facd*(1.-gl[k])*d1m1[index_mu][k])
                                     *( X_132_hhbar*X_121_hhbar
                                       +X_132_hhbar*X_121_hhbar)*d31[index_mu][l]
                                     +2.*(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *( X_132_hhbar_derv*X_121_hhbar
                                       +X_132_hhbar*X_121_hhbar_derv)*d31[index_mu][l]
                                     +2.*(facd*(1.-gl[k])*d1m1[index_mu][k])
                                     *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *( X_p022_hhbar*X_p022_hhbar*d22[index_mu][l]
                                       +X_242_hhbar*X_220_hhbar*d40[index_mu][l])
                                     +(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *(2.*X_p022_hhbar*X_p022_hhbar_derv*d22[index_mu][l]
                                       +( X_242_hhbar_derv*X_220_hhbar
                                         +X_242_hhbar*X_220_hhbar_derv)*d40[index_mu][l]));
                                    
                                    lensp_dl_derv += hl_P[l]*hl_P[l]*
                                    (2.*X_022_hh*X_022_hh_derv*d22[index_mu][l]+
                                     2.*(facd*(1.-2.*gl[k]+gl[k]*gl[k])
                                         *d1m1[index_mu][k])
                                     *X_132_hh*X_121_hh*d31[index_mu][l]
                                     +2.*( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                          +Cgl2_obs[index_mu])
                                     *( X_132_hh_derv*X_121_hh
                                       +X_132_hh*X_121_hh_derv)*d31[index_mu][l]
                                     +2.*(facd*(1.-2.*gl[k]+gl[k]*gl[k])
                                          *d1m1[index_mu][k])
                                     *( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *( X_p022_hh*X_p022_hh*d22[index_mu][l]
                                       +X_242_hh*X_220_hh*d40[index_mu][l])
                                     +( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *( 2.*X_p022_hh*X_p022_hh_derv*d22[index_mu][l]
                                       +( X_242_hh_derv*X_220_hh
                                         +X_242_hh*X_220_hh_derv)*d40[index_mu][l]));
                                }
                                
                                
                                /*------------ Cl_EE & CL_BB ------------*/
                                
                                resm_dl_derv = fac1*(cl_ee[l]-cl_bb[l]);
                                
                                lensm_dl_derv = 0.;
                                
                                lensm_dl_derv += hlbar_P[l]*hlbar_P[l]*
                                ( 2.*X_022*X_022_derv*d2m2[index_mu][l]
                                 +facd*d1m1[index_mu][k]
                                 *( X_121*X_121*d1m1[index_mu][l]
                                   +X_132*X_132*d3m3[index_mu][l])
                                 +Cgl2[index_mu]
                                 *( 2.*X_121*X_121_derv*d1m1[index_mu][l]
                                   +2.*X_132*X_132_derv*d3m3[index_mu][l])
                                 +facd*d1m1[index_mu][k]*Cgl2[index_mu]
                                 *( 2.*X_p022*X_p022*d2m2[index_mu][l]
                                   +X_220*X_220*d00[index_mu][l]
                                   +X_242*X_242*d4m4[index_mu][l])
                                 +Cgl2[index_mu]*Cgl2[index_mu]
                                 *( 2.*X_p022*X_p022_derv*d2m2[index_mu][l]
                                   +X_220*X_220_derv*d00[index_mu][l]
                                   +X_242*X_242_derv*d4m4[index_mu][l]));
                                
                                if(ple->derv_type == lensed)lensm_dl_derv /= hlbar_P[l]*hlbar_P[l];
                                
                                else if(ple->derv_type == delensed)
                                {
                                    lensm_dl_derv += 2.*hlbar_P[l]*hl_P[l]*
                                    (2.*X_022_hhbar*X_022_hhbar_derv*d2m2[index_mu][l]
                                     +(facd*(1.-gl[k])*d1m1[index_mu][k])
                                     *( X_121_hhbar*X_121_hhbar*d1m1[index_mu][l]
                                       +X_132_hhbar*X_132_hhbar*d3m3[index_mu][l])
                                     +(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *( 2.*X_121_hhbar*X_121_hhbar_derv*d1m1[index_mu][l]
                                       +2.*X_132_hhbar*X_132_hhbar_derv*d3m3[index_mu][l])
                                     +(facd*(1.-gl[k])*d1m1[index_mu][k])
                                     *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *( 2.*X_p022_hhbar*X_p022_hhbar*d2m2[index_mu][l]
                                       +X_220_hhbar*X_220_hhbar*d00[index_mu][l]
                                       +X_242_hhbar*X_242_hhbar*d4m4[index_mu][l])
                                     +(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *(Cgl2[index_mu]-Cgl2_cross[index_mu])
                                     *( 2.*X_p022_hhbar*X_p022_hhbar_derv*d2m2[index_mu][l]
                                       +X_220_hhbar*X_220_hhbar_derv*d00[index_mu][l]
                                       +X_242_hhbar*X_242_hhbar_derv*d4m4[index_mu][l]));
                                    
                                    lensm_dl_derv +=  hl_P[l]*hl_P[l]*
                                    (2.*X_022_hh*X_022_hh_derv*d2m2[index_mu][l]
                                     +(facd*(1.-2.*gl[k]+gl[k]*gl[k])
                                       *d1m1[index_mu][k])
                                     *(X_121_hh*X_121_hh*d1m1[index_mu][l]
                                       +X_132_hh*X_132_hh*d3m3[index_mu][l])
                                     +(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *( 2.*X_121_hh*X_121_hh_derv*d1m1[index_mu][l]
                                       +2.*X_132_hh*X_132_hh_derv*d3m3[index_mu][l])
                                     +(facd*(1.-2.*gl[k]+gl[k]*gl[k])
                                       *d1m1[index_mu][k])
                                     *(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *(2.*X_p022_hh*X_p022_hh*d2m2[index_mu][l]
                                       +X_220_hh*X_220_hh*d00[index_mu][l]
                                       +X_242_hh*X_242_hh*d4m4[index_mu][l])
                                     +( Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *(Cgl2[index_mu]-2.*Cgl2_cross[index_mu]
                                       +Cgl2_obs[index_mu])
                                     *(2.*X_p022_hh*X_p022_hh_derv*d2m2[index_mu][l]
                                       +X_220_hh*X_220_hh_derv*d00[index_mu][l]
                                       +X_242_hh*X_242_hh_derv*d4m4[index_mu][l]));
                                }
                                
                                if (ppr->accurate_lensing == _FALSE_) {
                                    lens_dl_derv  -= d00[index_mu][l];
                                    lensX_dl_derv -= d20[index_mu][l];
                                    lensm_dl_derv -= d22[index_mu][l];
                                    lensp_dl_derv -= d22[index_mu][l];
                                }
                                
                                res_dl_derv  *= lens_dl_derv;
                                resX_dl_derv *= lensX_dl_derv;
                                resm_dl_derv *= lensm_dl_derv;
                                resp_dl_derv *= lensp_dl_derv;
                                
                                ksi_dl_derv[index1_l][index_mu]  += res_dl_derv;
                                ksiX_dl_derv[index1_l][index_mu] += resX_dl_derv;
                                ksim_dl_derv[index1_l][index_mu] += resm_dl_derv;
                                ksip_dl_derv[index1_l][index_mu] += resp_dl_derv;
                                
                                
                            }
                        }
                    }
                }
            }

            if (ple->has_dl_tt == _TRUE_) {
                class_call(lensing_delensed_cl_tt(ksi_dl,hlbar,d00,w8,num_mu-1,ple),
                           ple->error_message,
                           ple->error_message);
                if (ppr->accurate_lensing == _FALSE_) {
                    class_call(lensing_addback_cl_dl_tt(ple,cl_tt),
                               ple->error_message,
                               ple->error_message);
                }
                // derivatives of the delensed spectra w.r.t. unlensed spectra
                /*
                if (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->delensed_wrt_unlensed == _TRUE_){
                    class_call(lensing_delensed_cl_tt_derv(ksi_dlu_derv,d00,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
                }
                */
                
                if(((ple->convergence_type == total && itr_index > 1
                      && cratio <= ple->convergence_criterion_itr) ||
                     (ple->convergence_type == every && itr_index > 1
                      && type_2_flag == _TRUE_) ||
                     ple->lens_rec_noise_type == external_rn ||
                     ple->has_itr_delensing == _FALSE_) &&
                    (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->delensed_wrt_unlensed == _TRUE_))
                {
                    class_call(lensing_delensed_cl_tt_derv_all(ksi_dlu_derv,d00,w8,num_mu-1,ple),
                           ple->error_message,
                           ple->error_message);
                }
                
                if(((ple->convergence_type == total && itr_index > 1
                     && cratio <= ple->convergence_criterion_itr) ||
                    (ple->convergence_type == every && itr_index > 1
                     && type_2_flag == _TRUE_) ||
                    ple->lens_rec_noise_type == external_rn || 
                    ple->has_itr_delensing == _FALSE_) &&
                   ple->calculate_pderivaties == _TRUE_)
                {
                    class_call(lensing_delensed_derv_cl_tt(ksi_dl_derv,hlbar,d00,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
                    if (ppr->accurate_lensing == _FALSE_) {
                        class_call(lensing_addback_derv_cl_dl_tt(ple,cl_tt),
                                   ple->error_message,
                                   ple->error_message);
                    }
                    

                    
                }
            }
            
            if (ple->has_dl_te == _TRUE_) {
                /*JRM: TE filtering
                 h filters in these functions do not do anything, and are therefore left unchanged */
                class_call(lensing_delensed_cl_te(ksiX_dl,hlbar,hlbar_P,d20,w8,num_mu-1,ple),
                           ple->error_message,
                           ple->error_message);
                if (ppr->accurate_lensing == _FALSE_) {
                    class_call(lensing_addback_cl_te(ple,cl_te),
                               ple->error_message,
                               ple->error_message);
                }
                /*
                if (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->delensed_wrt_unlensed == _TRUE_){
                    class_call(lensing_delensed_cl_te_derv(ksiX_dlu_derv,d20,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
                }
                */
                if(((ple->convergence_type == total && itr_index > 1
                      && cratio <= ple->convergence_criterion_itr) ||
                     (ple->convergence_type == every && itr_index > 1
                      && type_2_flag == _TRUE_) ||
                     ple->lens_rec_noise_type == external_rn ||
                     ple->has_itr_delensing == _FALSE_) &&
                    (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->delensed_wrt_unlensed == _TRUE_))
                {
                    class_call(lensing_delensed_cl_te_derv_all(ksiX_dlu_derv,d20,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
                }
                
                if(((ple->convergence_type == total && itr_index > 1
                     && cratio <= ple->convergence_criterion_itr) ||
                    (ple->convergence_type == every && itr_index > 1
                     && type_2_flag == _TRUE_) ||
                    ple->lens_rec_noise_type == external_rn || 
                    ple->has_itr_delensing == _FALSE_) &&
                   ple->calculate_pderivaties == _TRUE_)
                {
                    class_call(lensing_delensed_derv_cl_te(ksiX_dl_derv,hlbar,hlbar_P,d20,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
                    if (ppr->accurate_lensing == _FALSE_) {
                        class_call(lensing_addback_derv_cl_dl_te(ple,cl_te),
                                   ple->error_message,
                                   ple->error_message);
                    }
                    

                }
            }
            if (ple->has_dl_ee == _TRUE_ || ple->has_dl_bb == _TRUE_) {
                class_call(lensing_delensed_cl_ee_bb(ksip_dl,ksim_dl,hlbar_P,d22,d2m2,w8,num_mu-1,ple),
                           ple->error_message,
                           ple->error_message);
                if (ppr->accurate_lensing == _FALSE_) {
                    class_call(lensing_addback_cl_dl_ee_bb(ple,cl_ee,cl_bb),
                               ple->error_message,
                               ple->error_message);
                }
                
                /*
                if(ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->delensed_wrt_unlensed == _TRUE_){
//                     printf("above class_call lensing_delensed_cl_ee_bb_dervE \n");
                    class_call(lensing_delensed_cl_ee_bb_dervE(ksip_dlu_dervE,ksim_dlu_dervE,d22,d2m2,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
//                     printf("above class_call lensing_delensed_cl_ee_bb_dervB \n");
                    class_call(lensing_delensed_cl_ee_bb_dervB(ksip_dlu_dervB,ksim_dlu_dervB,d22,d2m2,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
                }
                */
                if(((ple->convergence_type == total && itr_index > 1
                      && cratio <= ple->convergence_criterion_itr) ||
                     (ple->convergence_type == every && itr_index > 1
                      && type_2_flag == _TRUE_) ||
                     ple->lens_rec_noise_type == external_rn ||
                     ple->has_itr_delensing == _FALSE_) &&
                    (ple->calculate_derviaties_wrt_unlensed == _TRUE_ && ple->delensed_wrt_unlensed == _TRUE_))
                {
                    class_call(lensing_delensed_cl_ee_bb_dervE_all(ksip_dlu_dervE,ksim_dlu_dervE,d22,d2m2,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
                    class_call(lensing_delensed_cl_ee_bb_dervB_all(ksip_dlu_dervB,ksim_dlu_dervB,d22,d2m2,w8,num_mu-1,ple),
                               ple->error_message,
                               ple->error_message);
                }
                
                 if(((ple->convergence_type == total && itr_index > 1
                      && cratio <= ple->convergence_criterion_itr) ||
                     (ple->convergence_type == every && itr_index > 1
                      && type_2_flag == _TRUE_) ||
                     ple->lens_rec_noise_type == external_rn ||
                     ple->has_itr_delensing == _FALSE_) &&
                    ple->calculate_pderivaties == _TRUE_)
                 {
                     class_call(lensing_delensed_derv_cl_ee_bb(ksip_dl_derv,ksim_dl_derv,hlbar_P,d22,d2m2,w8,num_mu-1,ple),
                                ple->error_message,
                                ple->error_message);
                     if (ppr->accurate_lensing == _FALSE_) {
                         class_call(lensing_addback_derv_cl_dl_ee_bb(ple,cl_ee,cl_bb),
                                    ple->error_message,
                                    ple->error_message);
                     }
                     

                 }
                
            }
            class_call(array_spline_table_lines(ple->l_dl,
                                                ple->dl_size,
                                                ple->cl_delens,
                                                ple->dlt_size,
                                                ple->ddcl_delens,
                                                _SPLINE_EST_DERIV_,
                                                ple->error_message),
                       ple->error_message,
                       ple->error_message);
            
            if(((ple->convergence_type == total && itr_index > 1
                 && cratio <= ple->convergence_criterion_itr) ||
                (ple->convergence_type == every && itr_index > 1
                 && type_2_flag == _TRUE_) ||
                ple->lens_rec_noise_type == external_rn || 
                ple->has_itr_delensing == _FALSE_) &&
               ple->calculate_pderivaties == _TRUE_)
            {
                if (ple->delensing_verbose > 0 && ple->derv_type == delensed)
                    printf("Calculating the 2nd-derivatives of the delensed spectra for sampling.\n");
                else if (ple->delensing_verbose > 0 && ple->derv_type == lensed)
                    printf("Calculating the 2nd-derivatives of the lensed spectra for sampling.\n");
                
                /*Given an m by n tabulated function ya[1..m][1..n], and tabulated independent variables
                 x2a[1..n], this routine constructs one-dimensional natural cubic splines of the rows of ya
                 and returns the second-derivatives in the array y2a[1..m][1..n].*/
                
                class_call(dlm_splie2(ple->l_dl,
                                      ple->l_dl,
                                      ple->cl_dl_tt_pderv,
                                      ple->dl_size-1,
                                      ple->dl_size-1,
                                      ple->ddcl_dl_tt_pderv),
                           ple->error_message,
                           ple->error_message);
                
                
                class_call(dlm_splie2(ple->l_dl,
                                      ple->l_dl,
                                      ple->cl_dl_te_pderv,
                                      ple->dl_size-1,
                                      ple->dl_size-1,
                                      ple->ddcl_dl_te_pderv),
                           ple->error_message,
                           ple->error_message);
                
                
                class_call(dlm_splie2(ple->l_dl,
                                      ple->l_dl,
                                      ple->cl_dl_ee_pderv,
                                      ple->dl_size-1,
                                      ple->dl_size-1,
                                      ple->ddcl_dl_ee_pderv),
                           ple->error_message,
                           ple->error_message);
                
                
                
                class_call(dlm_splie2(ple->l_dl,
                                      ple->l_dl,
                                      ple->cl_dl_bb_pderv,
                                      ple->dl_size-1,
                                      ple->dl_size-1,
                                      ple->ddcl_dl_bb_pderv),
                           ple->error_message,
                           ple->error_message);
                
                /*  
                //No need to spline when calculating at every l
                if(ple->calculate_derviaties_wrt_unlensed == _TRUE_){
                    if(ple->lensed_wrt_unlensed == _TRUE_){
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_lens_derv_TT_TT,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_lens_derv_TT_TT),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_lens_derv_TE_TE,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_lens_derv_TE_TE),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_lens_derv_EE_EE,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_lens_derv_EE_EE),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_lens_derv_EE_BB,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_lens_derv_EE_BB),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_lens_derv_BB_EE,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_lens_derv_BB_EE),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_lens_derv_BB_BB,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_lens_derv_BB_BB),
                                   ple->error_message,
                                   ple->error_message);
                        if (ple->delensing_verbose > 2)  printf("we are done with setting up the ddcls for lens d unlensed");
                    }
                    if(ple->delensed_wrt_unlensed == _TRUE_){
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_delens_derv_TT_TT,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_delens_derv_TT_TT),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_delens_derv_TE_TE,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_delens_derv_TE_TE),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_delens_derv_EE_EE,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_delens_derv_EE_EE),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_delens_derv_EE_BB,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_delens_derv_EE_BB),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_delens_derv_BB_EE,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_delens_derv_BB_EE),
                                   ple->error_message,
                                   ple->error_message);
                        
                        class_call(dlm_splie2(ple->l,
                                              ple->l,
                                              ple->cl_delens_derv_BB_BB,
                                              ple->l_size-1,
                                              ple->l_size-1,
                                              ple->ddcl_delens_derv_BB_BB),
                                   ple->error_message,
                                   ple->error_message);
                        if (ple->delensing_verbose > 2)  printf("we are done with setting up the ddcls for delens d unlensed");
                    }
                }
                */
                
            }
            for (l=2; l<=ple->l_delensed_max; l++) {
                
                ll = (double)l;
                
                class_call(delensing_cl_at_l(ple,
                                             ll,
                                             cl_delensed),
                           ple->error_message,
                           ple->error_message);
                
                if (ple->has_dl_tt == _TRUE_ &&
                    (ple->has_nl_eb_itr == _FALSE_  ||
                     ple->has_itr_delensing == _TRUE_)){
                    cl_delens_tt[l] = cl_delensed[ple->index_lt_dl_tt];
                    
                }
                if (ple->has_dl_te == _TRUE_ &&
                    (ple->has_nl_eb_itr == _FALSE_  ||
                     ple->has_itr_delensing == _TRUE_)){
                    cl_delens_te[l] = cl_delensed[ple->index_lt_dl_te];
                    
                }
                if (ple->has_dl_bb == _TRUE_ || ple->has_dl_ee == _TRUE_){
                    
                    if(ple->has_nl_eb_itr == _FALSE_ ||
                       ple->has_itr_delensing == _TRUE_) cl_delens_ee[l] = cl_delensed[ple->index_lt_dl_ee];
                    
                    cl_delens_bb[l] = cl_delensed[ple->index_lt_dl_bb];
                    
                }
                
            }
            
            if(ple->has_itr_delensing == _FALSE_) itr_index = ple->max_itr_steps + 1;
            
            ple->has_lens_noise_rcn = _FALSE_;
        }
    }
    /*---- end of the main delensing extension --------*/
    
    /** - Free lots of stuff **/
    free(buf_dxx);
    
    free(d00);
    free(d11);
    free(d1m1);
    free(d2m2);
    
    if (ple->has_te==_TRUE_) {
        free(d20);
        free(d3m1);
        free(d4m2);
    }
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        free(d22);
        free(d31);
        free(d3m3);
        
        if (ple->has_delensed_cls == _TRUE_){ /* DLM */
            free(d33);  /* DLM */
            free(d01);  /* DLM */
            free(d30);  /* DLM */
            free(d2m1); /* DLM */
            free(d21);  /* DLM */
            free(d32);  /* DLM */
            free(d3m2); /* DLM */
        }
        free(d40);
        free(d4m4);
    }
    
    if (ple->has_tt==_TRUE_)
        free(ksi);
    if (ple->has_te==_TRUE_)
        free(ksiX);
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        free(ksip);
        free(ksim);
    }
    
    free(Cgl);
    free(Cgl2);
    free(sigma2);
    
    free(mu);
    free(w8);
    
    free(cl_unlensed);
    free(cl_tt);
    if (ple->has_te==_TRUE_)
        free(cl_te);
    if (ple->has_ee==_TRUE_ || ple->has_bb==_TRUE_) {
        free(cl_ee);
        free(cl_bb);
    }
    
    if (ple->has_dl_tt==_TRUE_ && ple->has_delensed_cls == _TRUE_){ /* DLM */
        
        free(ksi_dl);  /* DLM */
        free(ksiX_dl); /* DLM */
        free(ksim_dl); /* DLM */
        free(ksip_dl); /* DLM */
        
        free(Cgl_obs);   /* DLM */
        free(Cgl_cross); /* DLM */
        free(Cgl2_obs);  /* DLM */
        free(Cgl2_cross);/* DLM */
        
        free(hl);   /* DLM */
        free(hlbar);/* DLM */
        free(gl);   /* DLM */
        
        free(hl_P);   /* DLM */
        free(hlbar_P);/* DLM */
        
        free(sigma2_hh);   /* DLM */
        free(sigma2_hhbar);/* DLM */
        
        free(cl_pp_obs);  /* DLM */
        free(cl_pp_cross);/* DLM */
    }
    
    if (ple->calculate_pderivaties ==_TRUE_ && ple->has_delensed_cls == _TRUE_){
        
        free(buf2_dxx);
        free(buf3_dxx);
        free(buf4_dxx);
        
        free(ksi_dl_derv);  /* DLM */
        free(ksiX_dl_derv); /* DLM */
        free(ksim_dl_derv); /* DLM */
        free(ksip_dl_derv); /* DLM */
    }
    
    if(ple->calculate_derviaties_wrt_unlensed == _TRUE_){
    
        if(ple->lensed_wrt_unlensed == _TRUE_){
            free(ksi_ln_derv);  /* DLM */
            free(ksiX_ln_derv); /* DLM */
            free(ksip_ln_dervE); /* DLM */
            free(ksip_ln_dervB); /* DLM */
            free(ksim_ln_dervE); /* DLM */
            free(ksim_ln_dervB); /* DLM */
        }
        if(ple->delensed_wrt_unlensed == _TRUE_){
    
            free(ksi_dlu_derv);  /* DLM */
            free(ksiX_dlu_derv); /* DLM */
            free(ksim_dlu_dervE); /* DLM */
            free(ksip_dlu_dervE); /* DLM */
            free(ksim_dlu_dervB); /* DLM */
            free(ksip_dlu_dervB); /* DLM */
        }
    }
    
    free(cl_pp);
    /** - Exit **/
    
    return _SUCCESS_;
    
}





/**
 * This routine frees all the memory space allocated by lensing_init().
 *
 * To be called at the end of each run, only when no further calls to
 * lensing_cl_at_l() are needed.
 *
 * @param ple Input: pointer to lensing structure (which fields must be freed)
 * @return the error status
 */

int lensing_free(
                 struct lensing * ple
                 ) {
    if (ple->has_lensed_cls == _TRUE_) {
        
        free(ple->l);
        free(ple->cl_lens);
        free(ple->ddcl_lens);
        free(ple->l_max_lt);
        
        if (ple->has_delensed_cls ==_TRUE_){ /* DLM */
            
            free(ple->l_unlensed);         /* DLM */
            free(ple->l_dl);         /* DLM */
            free(ple->cl_delens);   /* DLM */
            free(ple->ddcl_delens); /* DLM */
            
            free(ple->l_tn);  /* DLM */
            free(ple->pk_tn); /* DLM */
            free(ple->l_pn);  /* DLM */
            free(ple->pk_pn); /* DLM */
            
            free(ple->nl_rcn);   /* DLM */
            free(ple->ddnl_rcn); /* DLM */
            
            if(ple->temperature_noise_type == external_tn){    /* DLM */
                free(ple->command_for_temp_noise_spec);       /* DLM */
            }
            if(ple->polarization_noise_type == external_pn){   /* DLM */
                free(ple->command_for_polarization_noise_spec); /* DLM */
            }
            if(ple->lens_rec_noise_type == external_rn)       /* DLM */
            {
                free(ple->l_rcn_ext);  /* DLM */
                free(ple->pk_rcn_ext); /* DLM */
                free(ple->command_for_lens_recon_noise_spec); /* DLM */
            }
            
            if(ple->calculate_pderivaties==_TRUE_){
                
                free(ple->cl_dl_tt_pderv);
                free(ple->cl_dl_te_pderv);
                free(ple->cl_dl_ee_pderv);
                free(ple->cl_dl_bb_pderv);
                
                free(ple->ddcl_dl_tt_pderv);
                free(ple->ddcl_dl_te_pderv);
                free(ple->ddcl_dl_ee_pderv);
                free(ple->ddcl_dl_bb_pderv);
                
                if(ple->calculate_derviaties_wrt_unlensed == _TRUE_){
                    if(ple->lensed_wrt_unlensed == _TRUE_){
                        free(ple->cl_lens_derv_TT_TT);
                        free(ple->cl_lens_derv_TE_TE);
                        free(ple->cl_lens_derv_EE_EE);
                        free(ple->cl_lens_derv_EE_BB);
                        free(ple->cl_lens_derv_BB_EE);
                        free(ple->cl_lens_derv_BB_BB);
                        
                        free(ple->ddcl_lens_derv_TT_TT);
                        free(ple->ddcl_lens_derv_TE_TE);
                        free(ple->ddcl_lens_derv_EE_EE);
                        free(ple->ddcl_lens_derv_EE_BB);
                        free(ple->ddcl_lens_derv_BB_EE);
                        free(ple->ddcl_lens_derv_BB_BB);
                    }
                    if(ple->delensed_wrt_unlensed == _TRUE_){
                        free(ple->cl_delens_derv_TT_TT);
                        free(ple->cl_delens_derv_TE_TE);
                        free(ple->cl_delens_derv_EE_EE);
                        free(ple->cl_delens_derv_EE_BB);
                        free(ple->cl_delens_derv_BB_EE);
                        free(ple->cl_delens_derv_BB_BB);
                        
                        free(ple->ddcl_delens_derv_TT_TT);
                        free(ple->ddcl_delens_derv_TE_TE);
                        free(ple->ddcl_delens_derv_EE_EE);
                        free(ple->ddcl_delens_derv_EE_BB);
                        free(ple->ddcl_delens_derv_BB_EE);
                        free(ple->ddcl_delens_derv_BB_BB);
                    }
                }
            }
        }
    }
    return _SUCCESS_;
}

/**
 * This routine defines indices and allocates tables in the lensing structure
 *
 * @param ppr  Input: pointer to precision structure
 * @param phr  Input: pointer to harmonic structure
 * @param ple  Input/output: pointer to lensing structure
 * @return the error status
 */

int lensing_indices(
                    struct precision * ppr,
                    struct harmonic * phr,
                    struct lensing * ple
                    ){
    
    int index_l;
    
    double ** cl_md_ic; /* array with argument
                         cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */
    
    double ** cl_md;    /* array with argument
                         cl_md[index_md][index_ct] */
    
    int index_md;
    int index_lt;
    
    /* indices of all Cl types (lensed and unlensed) */
    
    if (phr->has_tt == _TRUE_) {
        ple->has_tt = _TRUE_;
        ple->index_lt_tt=phr->index_ct_tt;
    }
    else {
        ple->has_tt = _FALSE_;
    }
    
    if (phr->has_ee == _TRUE_) {
        ple->has_ee = _TRUE_;
        ple->index_lt_ee=phr->index_ct_ee;
    }
    else {
        ple->has_ee = _FALSE_;
    }
    
    if (phr->has_te == _TRUE_) {
        ple->has_te = _TRUE_;
        ple->index_lt_te=phr->index_ct_te;
    }
    else {
        ple->has_te = _FALSE_;
    }
    
    if (phr->has_bb == _TRUE_) {
        ple->has_bb = _TRUE_;
        ple->index_lt_bb=phr->index_ct_bb;
    }
    else {
        ple->has_bb = _FALSE_;
    }
    
    if (phr->has_pp == _TRUE_) {
        ple->has_pp = _TRUE_;
        ple->index_lt_pp=phr->index_ct_pp;
    }
    else {
        ple->has_pp = _FALSE_;
    }
    
    if (phr->has_tp == _TRUE_) {
        ple->has_tp = _TRUE_;
        ple->index_lt_tp=phr->index_ct_tp;
    }
    else {
        ple->has_tp = _FALSE_;
    }
    
    if (phr->has_dd == _TRUE_) {
        ple->has_dd = _TRUE_;
        ple->index_lt_dd=phr->index_ct_dd;
    }
    else {
        ple->has_dd = _FALSE_;
    }
    
    if (phr->has_td == _TRUE_) {
        ple->has_td = _TRUE_;
        ple->index_lt_td=phr->index_ct_td;
    }
    else {
        ple->has_td = _FALSE_;
    }
    
    if (phr->has_ll == _TRUE_) {
        ple->has_ll = _TRUE_;
        ple->index_lt_ll=phr->index_ct_ll;
    }
    else {
        ple->has_ll = _FALSE_;
    }
    
    if (phr->has_tl == _TRUE_) {
        ple->has_tl = _TRUE_;
        ple->index_lt_tl=phr->index_ct_tl;
    }
    else {
        ple->has_tl = _FALSE_;
    }
    
    ple->lt_size = phr->ct_size;
    
    /*-----------------------------------*/
    /*---NEEDS WORK BEFORE MAKES SENSE---*/
    /*-----------------------------------*/
    ple->dlt_size = 0; /* DLM */
    
    if (ple->has_delensed_cls == _TRUE_) { /* DLM */
        
        if (phr->has_tt == _TRUE_) { /* DLM */
            ple->has_dl_tt = _TRUE_; /* DLM */
            ple->index_lt_dl_tt=phr->index_ct_tt; /* DLM */
        }
        else {
            ple->has_dl_tt = _FALSE_; /* DLM */
        }
        
        if (phr->has_ee == _TRUE_) { /* DLM */
            ple->has_dl_ee = _TRUE_; /* DLM */
            ple->index_lt_dl_ee=phr->index_ct_ee; /* DLM */
        }
        else {
            ple->has_dl_ee = _FALSE_; /* DLM */
        }
        
        if (phr->has_te == _TRUE_) { /* DLM */
            ple->has_dl_te = _TRUE_; /* DLM */
            ple->index_lt_dl_te=phr->index_ct_te; /* DLM */
        }
        else {
            ple->has_dl_te = _FALSE_; /* DLM */
        }
        
        if (phr->has_bb == _TRUE_) { /* DLM */
            ple->has_dl_bb = _TRUE_; /* DLM */
            ple->index_lt_dl_bb=phr->index_ct_bb; /* DLM */
        }
        else {
            ple->has_dl_bb = _FALSE_; /* DLM */
        }
        
        if (phr->has_pp == _TRUE_) { /* DLM */
            ple->has_dl_pp = _TRUE_; /* DLM */
            ple->index_lt_dl_pp=phr->index_ct_pp; /* DLM */
        }
        else {
            ple->has_dl_pp = _FALSE_; /* DLM */
        }
        
        if (phr->has_tp == _TRUE_) { /* DLM */
            ple->has_dl_tp = _TRUE_; /* DLM */
            ple->index_lt_dl_tp=phr->index_ct_tp; /* DLM */
        }
        else {
            ple->has_dl_tp = _FALSE_; /* DLM */
        }
        
        if (phr->has_dd == _TRUE_) { /* DLM */
            ple->has_dl_dd = _TRUE_; /* DLM */
            ple->index_lt_dl_dd=phr->index_ct_dd; /* DLM */
        }
        else {
            ple->has_dl_dd = _FALSE_; /* DLM */
        }
        
        if (phr->has_td == _TRUE_) { /* DLM */
            ple->has_dl_td = _TRUE_; /* DLM */
            ple->index_lt_dl_td=phr->index_ct_td; /* DLM */
        }
        else {
            ple->has_dl_td = _FALSE_; /* DLM */
        }
        
        if (phr->has_ll == _TRUE_) { /* DLM */
            ple->has_dl_ll = _TRUE_; /* DLM */
            ple->index_lt_dl_ll=phr->index_ct_ll; /* DLM */
        }
        else {
            ple->has_dl_ll = _FALSE_; /* DLM */
        }
        
        if (phr->has_tl == _TRUE_) { /* DLM */
            ple->has_dl_tl = _TRUE_; /* DLM */
            ple->index_lt_dl_tl=phr->index_ct_tl; /* DLM */
        }
        else {
            ple->has_dl_tl = _FALSE_; /* DLM */
        }
        
        ple->dlt_size = phr->ct_size; /* DLM */
    }
    
    ple->nlt_size = 0; /* DLM */
    
    if (ple->lens_rec_noise_type == internal_rn) /* DLM */
    {
        
        if (ple->has_nl_all == _TRUE_) /* DLM */
        {
            /* all noise varriances and covariances
             used to generate the minimum variance. */
            ple->index_nl_minvar=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_tt=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_te=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_ee=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_bb=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_eb=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_tb=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_ttte=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_ttee=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_teee=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_tbeb=ple->nlt_size;
            ple->nlt_size++;
            
            if(ple->has_nl_all_itr == _TRUE_) /* DLM */ // CAN REMOVE NOT USED FOR ANYTHING
            {
                ple->has_ln_rcn_tt = _TRUE_;
                ple->has_ln_rcn_te = _TRUE_;
                ple->has_ln_rcn_ee = _TRUE_;
                ple->has_ln_rcn_bb = _TRUE_;
                ple->has_ln_rcn_eb = _TRUE_;
                ple->has_ln_rcn_tb = _TRUE_;
                ple->has_ln_rcn_ttte = _TRUE_;
                ple->has_ln_rcn_ttee = _TRUE_;
                ple->has_ln_rcn_teee = _TRUE_;
                ple->has_ln_rcn_tbeb = _TRUE_;
            }
            else if(ple->has_nl_diag_itr == _TRUE_) /* DLM */ // CAN REMOVE NOT USED FOR ANYTHING
            {
                ple->has_ln_rcn_tt = _TRUE_;
                ple->has_ln_rcn_te = _TRUE_;
                ple->has_ln_rcn_ee = _TRUE_;
                ple->has_ln_rcn_bb = _TRUE_;
                ple->has_ln_rcn_eb = _TRUE_;
                ple->has_ln_rcn_tb = _TRUE_;
            }
            else if(ple->has_nl_altr_itr == _TRUE_) /* DLM */ // CAN REMOVE NOT USED FOR ANYTHING
            {
                ple->has_ln_rcn_eb = _TRUE_;
                ple->has_ln_rcn_tb = _TRUE_;
            }
            
        }
        else if(ple->has_nl_diag == _TRUE_) /* DLM */
        {
            
            /* only varriances for the lensing
             reconstruction noise are used. */
            ple->index_nl_minvar=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_tt=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_te=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_ee=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_bb=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_eb=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_tb=ple->nlt_size;
            ple->nlt_size++;
            
            if(ple->has_nl_diag_itr == _TRUE_) /* DLM */
            {
                ple->has_ln_rcn_tt = _TRUE_;
                ple->has_ln_rcn_te = _TRUE_;
                ple->has_ln_rcn_ee = _TRUE_;
                ple->has_ln_rcn_bb = _TRUE_;
                ple->has_ln_rcn_eb = _TRUE_;
                ple->has_ln_rcn_tb = _TRUE_;
            }
            else if(ple->has_nl_altr_itr == _TRUE_) /* DLM */
            {
                ple->has_ln_rcn_eb = _TRUE_;
                ple->has_ln_rcn_tb = _TRUE_;
            }
        }
        else if(ple->has_nl_eb == _TRUE_) /* DLM */
        {
            
            /* only varriances for the lensing
             reconstruction noise are used. */
            ple->index_nl_minvar=ple->nlt_size;
            ple->nlt_size++;
            ple->index_nl_eb=ple->nlt_size;
            ple->nlt_size++;
        }
        
        // <<<<< COMPLETE WITH OTHER(?) OPTIONS
    }
    else{
        ple->index_nl_minvar=ple->nlt_size;
        ple->nlt_size++;
    }
    
    /* number of multipoles */
    
    ple->l_unlensed_max = phr->l_max_tot;
    
    ple->l_lensed_max = ple->l_unlensed_max - ppr->delta_l_max;
    
    for (index_l=0; (index_l < phr->l_size_max) && (phr->l[index_l] <= ple->l_lensed_max); index_l++);
    
    if (index_l < phr->l_size_max) index_l++; /* one more point in order to be able to interpolate till ple->l_lensed_max */
    
    ple->l_size = index_l; /* DLM */
    
    class_alloc(ple->l,ple->l_size*sizeof(double),ple->error_message);
    
    for (index_l=0; index_l < ple->l_size; index_l++) {
        
        ple->l[index_l] = phr->l[index_l];
        
    }
    
    class_alloc(ple->l_unlensed,(ple->l_unlensed_max+1)*sizeof(double),ple->error_message);
    
    for (index_l=0; index_l < ple->l_unlensed_max+1; index_l++) {
        
        ple->l_unlensed[index_l] = (double)index_l;
        
    }
    
    
    ple->l_delensed_max = ple->l_lensed_max - ppr->delta_dl_max; /* DLM */
    
    for (index_l=0; (index_l < phr->l_size_max) && (phr->l[index_l] <= ple->l_delensed_max); index_l++);  /* DLM */
    
    if (index_l < phr->l_size_max) index_l++; /* DLM */
    
    ple->dl_size = index_l; /* DLM */
    
    class_alloc(ple->l_dl,ple->dl_size*sizeof(double),ple->error_message);  /* DLM */
    
    for (index_l=0; index_l < ple->dl_size; index_l++) {  /* DLM */
        
        ple->l_dl[index_l] = phr->l[index_l]; /* DLM */
        
    }
    
    /* allocate table where results will be stored */
    
    class_alloc(ple->cl_lens,
                ple->l_size*ple->lt_size*sizeof(double),
                ple->error_message);
    
    class_alloc(ple->ddcl_lens,
                ple->l_size*ple->lt_size*sizeof(double),
                ple->error_message);
    
    
    if(ple->has_delensed_cls ==_TRUE_){ /* DLM */
        
        class_alloc(ple->nl_rcn,
                    ple->dl_size*ple->nlt_size*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(ple->ddnl_rcn,
                    ple->dl_size*ple->nlt_size*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(ple->cl_delens,
                    ple->dl_size*ple->dlt_size*sizeof(double),
                    ple->error_message); /* DLM */
        
        class_alloc(ple->ddcl_delens,
                    ple->dl_size*ple->dlt_size*sizeof(double),
                    ple->error_message); /* DLM */
        
    }
    
    /* fill with unlensed cls */
    
    class_alloc(cl_md_ic,
                phr->md_size*sizeof(double *),
                ple->error_message);
    
    class_alloc(cl_md,
                phr->md_size*sizeof(double *),
                ple->error_message);
    
    for (index_md = 0; index_md < phr->md_size; index_md++) {
        
        if (phr->md_size > 1)
            
            class_alloc(cl_md[index_md],
                        phr->ct_size*sizeof(double),
                        ple->error_message);
        
        if (phr->ic_size[index_md] > 1)
            
            class_alloc(cl_md_ic[index_md],
                        phr->ic_ic_size[index_md]*phr->ct_size*sizeof(double),
                        ple->error_message);
    }
    
    for (index_l=0; index_l<ple->l_size; index_l++) {
        
        class_call(harmonic_cl_at_l(phr,ple->l[index_l],&(ple->cl_lens[index_l*ple->lt_size]),cl_md,cl_md_ic),
                   phr->error_message,
                   ple->error_message);
        
    }
    
    /*----------------SH--------------*/
    /*---- probabily unnecessary -----*/
    if (ple->has_delensed_cls == _TRUE_) {
        
        for (index_l=0; index_l<ple->dl_size; index_l++) {
            
            
            if( ple->has_dl_tt == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_tt] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_tt];
            
            if( ple->has_dl_pp == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_pp] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_pp];
            
            if( ple->has_dl_tp == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_tp] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_tp];
            
            if( ple->has_dl_te == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_te] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_te];
            
            if( ple->has_dl_ee == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_ee] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_ee];
            
            if( ple->has_dl_bb == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_bb] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_bb];
            
            if( ple->has_dl_dd == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_dd] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_dd];
            
            if( ple->has_dl_td == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_td] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_td];
            
            if( ple->has_dl_ll == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_ll] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_ll];
            
            if( ple->has_dl_tl == _TRUE_ ) ple->cl_delens[index_l*ple->dlt_size+
                                                          ple->index_lt_dl_tl] = ple->cl_lens[index_l*ple->lt_size+
                                                                                              ple->index_lt_tl];
        }
    }
    /*--------------------------------*/
    
    for (index_md = 0; index_md < phr->md_size; index_md++) {
        
        if (phr->md_size > 1)
            free(cl_md[index_md]);
        
        if (phr->ic_size[index_md] > 1)
            free(cl_md_ic[index_md]);
        
    }
    
    free(cl_md_ic);
    free(cl_md);
    
    /* we want to output Cl_lensed up to the same l_max as Cl_unlensed
     (even if a number delta_l_max of extra values of l have been used
     internally for more accurate results). Notable exception to the
     above rule: ClBB_lensed(scalars) must be outputed at least up to the same l_max as
     ClEE_unlensed(scalars) (since ClBB_unlensed is null for scalars)
     */
    
    class_alloc(ple->l_max_lt,ple->lt_size*sizeof(double),ple->error_message);
    for (index_lt = 0; index_lt < ple->lt_size; index_lt++) {
        ple->l_max_lt[index_lt]=0.;
        for (index_md = 0; index_md < phr->md_size; index_md++) {
            ple->l_max_lt[index_lt]=MAX(ple->l_max_lt[index_lt],phr->l_max_ct[index_md][index_lt]);
            
            if ((ple->has_bb == _TRUE_) && (ple->has_ee == _TRUE_) && (index_lt == ple->index_lt_bb)) {
                ple->l_max_lt[index_lt]=MAX(ple->l_max_lt[index_lt],phr->l_max_ct[index_md][ple->index_lt_ee]);
            }
            
        }
    }
    
    return _SUCCESS_;
    
}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksi  Input: Lensed correlation function (ksi[index_mu])
 * @param d00  Input: Legendre polynomials (\f$ d^l_{00}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */


int lensing_lensed_cl_tt(
                         double *ksi,
                         double **d00,
                         double *w8,
                         int nmu,
                         struct lensing * ple
                         ) {
    
    double cle;
    int imu;
    int index_l;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l,cle)                     \
schedule (static)
    
    for(index_l=0; index_l<ple->l_size; index_l++){
        cle=0;
        for (imu=0;imu<nmu;imu++) {
            cle += ksi[imu]*d00[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
        }
        ple->cl_lens[index_l*ple->lt_size+ple->index_lt_tt]=cle*2.0*_PI_;
    }
    
    return _SUCCESS_;
}

int lensing_lensed_cl_tt_derv(
                              double **ksi,
                              double **d00,
                              double *w8,
                              int nmu,
                              struct lensing * ple
                              ) {
    
    double cle;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,cle)         \
schedule (static)
    
    for(index_l_1=0; index_l_1<ple->l_size; index_l_1++){
        for(index_l_2=0; index_l_2<ple->l_size; index_l_2++){
            cle=0;
            for (imu=0;imu<nmu;imu++) {
                cle += ksi[(int)ple->l[index_l_2]][imu]*d00[imu][(int)ple->l[index_l_1]]*w8[imu]; /* loop could be optimized */
            }
            ple->cl_lens_derv_TT_TT[index_l_2][index_l_1]=cle*2.0*_PI_;
        }
    }
    return _SUCCESS_;
}

int lensing_lensed_cl_tt_derv_all(
                              double **ksi,
                              double **d00,
                              double *w8,
                              int nmu,
                              struct lensing * ple
                              ) {
    
    double cle;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,cle)         \
schedule (static)
    
    for(index_l_1=2; index_l_1<ple->l_lensed_max+1; index_l_1++){
        for(index_l_2=2; index_l_2<ple->l_lensed_max+1; index_l_2++){
            cle=0;
            for (imu=0;imu<nmu;imu++) {
                cle += ksi[index_l_2][imu]*d00[imu][index_l_1]*w8[imu]; /* loop could be optimized */
            }
            // dCl^TT/dCL^{TT,unlensed}
            // index_l_1 = l, index_l_2 = L
            // Index ordering here is chosen for contiguous memory blocks during output
            ple->cl_lens_derv_TT_TT[index_l_1][index_l_2]=cle*2.0*_PI_;
        }
    }
    return _SUCCESS_;
}

  for (index_l=0; index_l<ple->l_size; index_l++){
    cle=0;
    for (imu=0;imu<nmu;imu++) {
      cle += ksi[imu]*d00[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    return _SUCCESS_;
}

int lensing_delensed_cl_tt_derv_all(
                              double **ksi,
                              double **d00,
                              double *w8,
                              int nmu,
                              struct lensing * ple
                              ) {
    
    double cle;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,cle)         \
schedule (static)
    
    for(index_l_1=2; index_l_1<ple->l_lensed_max+1; index_l_1++){
        for(index_l_2=2; index_l_2<ple->l_lensed_max+1; index_l_2++){
            cle=0;
            for (imu=0;imu<nmu;imu++) {
                cle += ksi[index_l_2][imu]*d00[imu][index_l_1]*w8[imu]; /* loop could be optimized */
            }
            ple->cl_delens_derv_TT_TT[index_l_1][index_l_2]=cle*2.0*_PI_;
        }
    }
    return _SUCCESS_;
}

/**
 * This routine adds back the unlensed \f$ cl_{tt}\f$ power spectrum
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param ple   Input/output: Pointer to the lensing structure
 * @param cl_tt Input: Array of unlensed power spectrum
 * @return the error status
 */

int lensing_addback_cl_tt(
                          struct lensing * ple,
                          double *cl_tt) {
    int index_l, l;
    
    for (index_l=0; index_l<ple->l_size; index_l++) {
        l = (int)ple->l[index_l];
        ple->cl_lens[index_l*ple->lt_size+ple->index_lt_tt] += cl_tt[l];
    }
    return _SUCCESS_;
    
}
// think about addback functions for lensing dervs
//int lensing_addback_cl_tt_derv(
//                          struct lensing * ple,
//                          double *cl_tt) {
//  int index_l, l;
//
//  for (index_l=0; index_l<ple->l_size; index_l++) {
//    l = (int)ple->l[index_l];
//    ple->cl_lens_derv[index_l*ple->lt_size+ple->index_lt_tt] += 1.;
//  }
//  return _SUCCESS_;
//
//}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksiX Input: Lensed correlation function (ksiX[index_mu])
 * @param d20  Input: Wigner d-function (\f$ d^l_{20}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */


int lensing_lensed_cl_te(
                         double *ksiX,
                         double **d20,
                         double *w8,
                         int nmu,
                         struct lensing * ple
                         ) {
    
    double clte;
    int imu;
    int index_l;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l,clte)                    \
schedule (static)
    
    for(index_l=0; index_l < ple->l_size; index_l++){
        clte=0;
        for (imu=0;imu<nmu;imu++) {
            clte += ksiX[imu]*d20[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
        }
        ple->cl_lens[index_l*ple->lt_size+ple->index_lt_te]=clte*2.0*_PI_;
    }
    
    return _SUCCESS_;
}

int lensing_lensed_cl_te_derv(
                              double **ksiX,
                              double **d20,
                              double *w8,
                              int nmu,
                              struct lensing * ple
                              ) {
    
    double clte;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clte)         \
schedule (static)
    
    for(index_l_1=0; index_l_1<ple->l_size; index_l_1++){
        for(index_l_2=0; index_l_2<ple->l_size; index_l_2++){
            clte=0;
            for (imu=0;imu<nmu;imu++) {
                clte += ksiX[(int)ple->l[index_l_2]][imu]*d20[imu][(int)ple->l[index_l_1]]*w8[imu]; /* loop could be optimized */
            }
            ple->cl_lens_derv_TE_TE[index_l_2][index_l_1]=clte*2.0*_PI_;
        }
    }
    return _SUCCESS_;
}

int lensing_lensed_cl_te_derv_all(
                              double **ksiX,
                              double **d20,
                              double *w8,
                              int nmu,
                              struct lensing * ple
                              ) {
    
    double clte;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clte)         \
schedule (static)
    
    for(index_l_1=2; index_l_1<ple->l_lensed_max+1; index_l_1++){
        for(index_l_2=2; index_l_2<ple->l_lensed_max+1; index_l_2++){
            clte=0;
            for (imu=0;imu<nmu;imu++) {
                clte += ksiX[index_l_2][imu]*d20[imu][index_l_1]*w8[imu]; /* loop could be optimized */
            }
            ple->cl_lens_derv_TE_TE[index_l_1][index_l_2]=clte*2.0*_PI_;
        }
    }
    return _SUCCESS_;
}

  for (index_l=0; index_l < ple->l_size; index_l++){
    clte=0;
    for (imu=0;imu<nmu;imu++) {
      clte += ksiX[imu]*d20[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    return _SUCCESS_;
}

int lensing_delensed_cl_te_derv_all(
                              double **ksiX,
                              double **d20,
                              double *w8,
                              int nmu,
                              struct lensing * ple
                              ) {
    
    double clte;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clte)         \
schedule (static)
    
    for(index_l_1=2; index_l_1<ple->l_lensed_max+1; index_l_1++){
        for(index_l_2=2; index_l_2<ple->l_lensed_max+1; index_l_2++){
            clte=0;
            for (imu=0;imu<nmu;imu++) {
                clte += ksiX[index_l_2][imu]*d20[imu][index_l_1]*w8[imu]; /* loop could be optimized */
            }
            ple->cl_delens_derv_TE_TE[index_l_1][index_l_2]=clte*2.0*_PI_;
        }
    }
    return _SUCCESS_;
}


/**
 * This routine adds back the unlensed \f$ cl_{te}\f$ power spectrum
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param ple   Input/output: Pointer to the lensing structure
 * @param cl_te Input: Array of unlensed power spectrum
 * @return the error status
 */

int lensing_addback_cl_te(
                          struct lensing * ple,
                          double *cl_te) {
    int index_l, l;
    
    for (index_l=0; index_l<ple->l_size; index_l++) {
        l = (int)ple->l[index_l];
        ple->cl_lens[index_l*ple->lt_size+ple->index_lt_te] += cl_te[l];
    }
    return _SUCCESS_;
    
}

//int lensing_addback_cl_te_derv(
//                          struct lensing * ple,
//                          double *cl_te) {
//  int index_l, l;
//
//  for (index_l=0; index_l<ple->l_size; index_l++) {
//    l = (int)ple->l[index_l];
//    ple->cl_lens_derv[index_l*ple->lt_size+ple->index_lt_te] += 1.;
//  }
//  return _SUCCESS_;
//
//}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksip Input: Lensed correlation function (ksi+[index_mu])
 * @param ksim Input: Lensed correlation function (ksi-[index_mu])
 * @param d22  Input: Wigner d-function (\f$ d^l_{22}\f$[l][index_mu])
 * @param d2m2 Input: Wigner d-function (\f$ d^l_{2-2}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */


int lensing_lensed_cl_ee_bb(
                            double *ksip,
                            double *ksim,
                            double **d22,
                            double **d2m2,
                            double *w8,
                            int nmu,
                            struct lensing * ple
                            ) {
    
    double clp, clm;
    int imu;
    int index_l;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l,clp,clm)                 \
schedule (static)
    
    for(index_l=0; index_l < ple->l_size; index_l++){
        clp=0; clm=0;
        for (imu=0;imu<nmu;imu++) {
            clp += ksip[imu]*d22[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
            clm += ksim[imu]*d2m2[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
        }
        ple->cl_lens[index_l*ple->lt_size+ple->index_lt_ee]=(clp+clm)*_PI_;
        ple->cl_lens[index_l*ple->lt_size+ple->index_lt_bb]=(clp-clm)*_PI_;
    }
    
    return _SUCCESS_;
}

int lensing_lensed_cl_ee_bb_dervE(
                                  double **ksip,
                                  double **ksim,
                                  double **d22,
                                  double **d2m2,
                                  double *w8,
                                  int nmu,
                                  struct lensing * ple
                                  ) {
    
    double clp, clm;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clp,clm)                 \
schedule (static)
    
    for(index_l_1=0; index_l_1<ple->dl_size; index_l_1++){
        for(index_l_2=0; index_l_2<ple->dl_size; index_l_2++){
            clp=0; clm=0;
            for (imu=0;imu<nmu;imu++) {
                clp += ksip[(int)ple->l_dl[index_l_2]][imu]*d22[imu][(int)ple->l_dl[index_l_1]]*w8[imu]; 
                clm += ksim[(int)ple->l_dl[index_l_2]][imu]*d2m2[imu][(int)ple->l_dl[index_l_1]]*w8[imu];
            }
            ple->cl_lens_derv_EE_EE[index_l_2][index_l_1]=(clp+clm)*_PI_;
            ple->cl_lens_derv_EE_BB[index_l_2][index_l_1]=(clp-clm)*_PI_;
        }
    }
    return _SUCCESS_;
}

int lensing_lensed_cl_ee_bb_dervE_all(
                                  double **ksip,
                                  double **ksim,
                                  double **d22,
                                  double **d2m2,
                                  double *w8,
                                  int nmu,
                                  struct lensing * ple
                                  ) {
    
    double clp, clm;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clp,clm)                 \
schedule (static)
    
    for(index_l_1=2; index_l_1<ple->l_lensed_max+1; index_l_1++){
        for(index_l_2=2; index_l_2<ple->l_lensed_max+1; index_l_2++){
            clp=0; clm=0;
            for (imu=0;imu<nmu;imu++) {
                clp += ksip[index_l_2][imu]*d22[imu][index_l_1]*w8[imu]; 
                clm += ksim[index_l_2][imu]*d2m2[imu][index_l_1]*w8[imu];
            }
            ple->cl_lens_derv_EE_EE[index_l_1][index_l_2]=(clp+clm)*_PI_;
            ple->cl_lens_derv_EE_BB[index_l_1][index_l_2]=(clp-clm)*_PI_;
        }
    }
    return _SUCCESS_;
}

  for (index_l=0; index_l < ple->l_size; index_l++){
    clp=0; clm=0;
    for (imu=0;imu<nmu;imu++) {
      clp += ksip[imu]*d22[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
      clm += ksim[imu]*d2m2[imu][(int)ple->l[index_l]]*w8[imu]; /* loop could be optimized */
    }
    return _SUCCESS_;
}

int lensing_delensed_cl_ee_bb_dervE_all(
                                  double **ksip,
                                  double **ksim,
                                  double **d22,
                                  double **d2m2,
                                  double *w8,
                                  int nmu,
                                  struct lensing * ple
                                  ) {
    
    double clp, clm;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clp,clm)                 \
schedule (static)
    
    for(index_l_1=2; index_l_1<ple->l_lensed_max+1; index_l_1++){
        for(index_l_2=2; index_l_2<ple->l_lensed_max+1; index_l_2++){
            clp=0; clm=0;
            for (imu=0;imu<nmu;imu++) {
                clp += ksip[index_l_2][imu]*d22[imu][index_l_1]*w8[imu]; 
                clm += ksim[index_l_2][imu]*d2m2[imu][index_l_1]*w8[imu];
            }
            ple->cl_delens_derv_EE_EE[index_l_1][index_l_2]=(clp+clm)*_PI_;
            ple->cl_delens_derv_EE_BB[index_l_1][index_l_2]=(clp-clm)*_PI_;
        }
    }
    return _SUCCESS_;
}

int lensing_lensed_cl_ee_bb_dervB(
                                  double **ksip,
                                  double **ksim,
                                  double **d22,
                                  double **d2m2,
                                  double *w8,
                                  int nmu,
                                  struct lensing * ple
                                  ) {
    
    double clp, clm;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clp,clm)                 \
schedule (static)
    
    for(index_l_1=0; index_l_1<ple->dl_size; index_l_1++){
        for(index_l_2=0; index_l_2<ple->dl_size; index_l_2++){
            clp=0; clm=0;
            for (imu=0;imu<nmu;imu++) {
                clp += ksip[(int)ple->l_dl[index_l_2]][imu]*d22[imu][(int)ple->l_dl[index_l_1]]*w8[imu];
                clm += ksim[(int)ple->l_dl[index_l_2]][imu]*d2m2[imu][(int)ple->l_dl[index_l_1]]*w8[imu]; 
            }
            ple->cl_lens_derv_BB_EE[index_l_2][index_l_1]=(clp+clm)*_PI_;
            ple->cl_lens_derv_BB_BB[index_l_2][index_l_1]=(clp-clm)*_PI_;
        }
    }
    return _SUCCESS_;
}

int lensing_lensed_cl_ee_bb_dervB_all(
                                  double **ksip,
                                  double **ksim,
                                  double **d22,
                                  double **d2m2,
                                  double *w8,
                                  int nmu,
                                  struct lensing * ple
                                  ) {
    
    double clp, clm;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clp,clm)                 \
schedule (static)
    
    for(index_l_1=2; index_l_1<ple->l_lensed_max+1; index_l_1++){
        for(index_l_2=2; index_l_2<ple->l_lensed_max+1; index_l_2++){
            clp=0; clm=0;
            for (imu=0;imu<nmu;imu++) {
                clp += ksip[index_l_2][imu]*d22[imu][index_l_1]*w8[imu];
                clm += ksim[index_l_2][imu]*d2m2[imu][index_l_1]*w8[imu]; 
            }
            ple->cl_lens_derv_BB_EE[index_l_1][index_l_2]=(clp+clm)*_PI_;
            ple->cl_lens_derv_BB_BB[index_l_1][index_l_2]=(clp-clm)*_PI_;
        }
    }
    return _SUCCESS_;
}

int lensing_delensed_cl_ee_bb_dervB_all(
                                  double **ksip,
                                  double **ksim,
                                  double **d22,
                                  double **d2m2,
                                  double *w8,
                                  int nmu,
                                  struct lensing * ple
                                  ) {
    
    double clp, clm;
    int imu;
    int index_l_1;
    int index_l_2;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clp,clm)                 \
schedule (static)
    
    for(index_l_1=2; index_l_1<ple->l_lensed_max+1; index_l_1++){
        for(index_l_2=2; index_l_2<ple->l_lensed_max+1; index_l_2++){
            clp=0; clm=0;
            for (imu=0;imu<nmu;imu++) {
                clp += ksip[index_l_2][imu]*d22[imu][index_l_1]*w8[imu];
                clm += ksim[index_l_2][imu]*d2m2[imu][index_l_1]*w8[imu]; 
            }
            ple->cl_delens_derv_BB_EE[index_l_1][index_l_2]=(clp+clm)*_PI_;
            ple->cl_delens_derv_BB_BB[index_l_1][index_l_2]=(clp-clm)*_PI_;
        }
    }
    return _SUCCESS_;
}



int lensing_delensed_cl_ee_bb_dervB(
                                    double **ksip,
                                    double **ksim,
                                    double **d22,
                                    double **d2m2,
                                    double *w8,
                                    int nmu,
                                    struct lensing * ple
                                    ) {
    
    double clp, clm;
    int imu;
    int index_l_1;
    int index_l_2;
    
//     printf("we are inside lensing_delensed_cl_ee_bb_dervB\n");    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l_1,index_l_2,clp,clm)                 \
schedule (static)

    for(index_l_1=0; index_l_1<ple->dl_size; index_l_1++){
        for(index_l_2=0; index_l_2<ple->dl_size; index_l_2++){
            clp=0; clm=0;
            for (imu=0;imu<nmu;imu++) {
                if(imu==0){
                }
                clp += ksip[(int)ple->l_dl[index_l_2]][imu]*d22[imu][(int)ple->l_dl[index_l_1]]*w8[imu]; /* loop could be optimized */
                clm += ksim[(int)ple->l_dl[index_l_2]][imu]*d2m2[imu][(int)ple->l_dl[index_l_1]]*w8[imu]; /* loop could be optimized */
            }
            ple->cl_delens_derv_BB_EE[index_l_2][index_l_1]=(clp+clm)*_PI_;
            ple->cl_delens_derv_BB_BB[index_l_2][index_l_1]=(clp-clm)*_PI_;
        }
    }
    return _SUCCESS_;
}

/**
 * This routine adds back the unlensed \f$ cl_{ee}\f$, \f$ cl_{bb}\f$ power spectra
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param ple   Input/output: Pointer to the lensing structure
 * @param cl_ee Input: Array of unlensed power spectrum
 * @param cl_bb Input: Array of unlensed power spectrum
 * @return the error status
 */

int lensing_addback_cl_ee_bb(
                             struct lensing * ple,
                             double * cl_ee,
                             double * cl_bb) {
    
    int index_l, l;
    
    for (index_l=0; index_l<ple->l_size; index_l++) {
        l = (int)ple->l[index_l];
        ple->cl_lens[index_l*ple->lt_size+ple->index_lt_ee] += cl_ee[l];
        ple->cl_lens[index_l*ple->lt_size+ple->index_lt_bb] += cl_bb[l];
    }
    return _SUCCESS_;
    
}

/**
 * DLM: This routine computes the delensed power spectra by Gaussian quadrature
 *
 * @param ksi_dl  Input: Delensed correlation function (ksi_dl[index_mu]) omitting lensed cl_tt term.
 * @param hl   Input: -----------
 * @param d00  Input: Legendre polynomials (\f$ d^l_{00}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */
int lensing_delensed_cl_tt(
                           double *ksi_dl,
                           double *hla,
                           double **d00,
                           double *w8,
                           int nmu,
                           struct lensing * ple
                           ) {
    
    double cdl,  ll; /* lensed part of the correleation function */
    int imu, index_l, l;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l,cdl) \
schedule (static)
    
    for(index_l=0; index_l<ple->dl_size; index_l++){
        
        cdl=0.;
        
        for (imu=0;imu<nmu;imu++) {
            
            cdl += ksi_dl[imu]*d00[imu][(int)ple->l_dl[index_l]]*w8[imu]; /* loop could be optimized */
            
        }
        
        cdl *= 2.0*_PI_; /* calculated the delensed part of the spectra*/
        
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_tt]=cdl;
        
    }
    
    return _SUCCESS_;
}

int lensing_addback_cl_dl_tt(
                             struct lensing * ple,
                             double *cl_tt) {
    int index_l, l;
    
    for (index_l=0; index_l<ple->l_size; index_l++) {
        l = (int)ple->l_dl[index_l];
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_tt] += cl_tt[l];
    }
    return _SUCCESS_;
    
}

int lensing_delensed_cl_te(
                           double *ksiX_dl,
                           double *hla,
                           double *hla_P,
                           double **d20,
                           double *w8,
                           int nmu,
                           struct lensing * ple
                           ) {
    
    double clte_dl;
    int imu;
    int index_l;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l,clte_dl)                    \
schedule (static)
    
    for(index_l=0; index_l < ple->dl_size; index_l++){
        
        clte_dl=0;
        
        for (imu=0;imu<nmu;imu++) {
            
            clte_dl += ksiX_dl[imu]*d20[imu][(int)ple->l_dl[index_l]]*w8[imu]; /* loop could be optimized */
            
        }
        
        clte_dl *= 2.0*_PI_;
        
        
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_te]=clte_dl;
        
    }
    
    return _SUCCESS_;
}

int lensing_addback_cl_dl_te(
                             struct lensing * ple,
                             double *cl_te) {
    int index_l, l;
    
    for (index_l=0; index_l<ple->l_size; index_l++) {
        l = (int)ple->l_dl[index_l];
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_te] += cl_te[l];
    }
    return _SUCCESS_;
    
}

int lensing_delensed_cl_ee_bb(
                              double *ksip_dl,
                              double *ksim_dl,
                              double *hla,
                              double **d22,
                              double **d2m2,
                              double *w8,
                              int nmu,
                              struct lensing * ple
                              ) {
    
    double clp_dl, clm_dl, cl_ee, cl_bb;
    int imu;
    int index_l;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l,clp_dl,clm_dl,cl_ee,cl_bb) \
schedule (static)
    
    for(index_l=0; index_l < ple->dl_size; index_l++){
        
        clp_dl=0; clm_dl=0;
        
        for (imu=0;imu<nmu;imu++) {
            
            clp_dl += ksip_dl[imu]*d22[imu][(int)ple->l_dl[index_l]]*w8[imu]; /* loop could be optimized */
            clm_dl += ksim_dl[imu]*d2m2[imu][(int)ple->l_dl[index_l]]*w8[imu]; /* loop could be optimized */
            
        }
        
        cl_ee = (clp_dl+clm_dl)*_PI_;
        cl_bb = (clp_dl-clm_dl)*_PI_;
        
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_ee]=cl_ee;
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_bb]=cl_bb;
        
        
    }
    
    return _SUCCESS_;
}

int lensing_addback_cl_dl_ee_bb(
                                struct lensing * ple,
                                double * cl_ee,
                                double * cl_bb) {
    
    int index_l, l;
    
    for (index_l=0; index_l<ple->dl_size; index_l++) {
        l = (int)ple->l_dl[index_l];
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_ee] += cl_ee[l];
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_bb] += cl_bb[l];
    }
    return _SUCCESS_;
}


/**
 * DLM: This routine computes the delensed power spectra by Gaussian quadrature
 *
 * @param ksi_dl  Input: Delensed correlation function (ksi_dl[index_mu]) omitting lensed cl_tt term.
 * @param hl   Input: -----------
 * @param d00  Input: Legendre polynomials (\f$ d^l_{00}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */
int lensing_delensed_derv_cl_tt(
                                double **ksi_dl_derv,
                                double *hla,
                                double **d00,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                ) {
    
    double cdl,  ll; /* lensed part of the correleation function */
    int imu, index2_l,index2_k,l;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index2_l,index2_k,cdl) \
schedule (static)
    
    for(index2_k=0; index2_k<ple->dl_size; index2_k++){
        for(index2_l=0; index2_l<ple->dl_size; index2_l++){
            
            cdl=0.;
            
            for (imu=0;imu<nmu;imu++) {
                
                cdl += ksi_dl_derv[index2_k][imu]*d00[imu][(int)ple->l_dl[index2_l]]*w8[imu]; /* loop could be optimized */
                
            }
            
            cdl *= 2.0*_PI_; /* calculated the delensed part of the spectra*/
            
            ple->cl_dl_tt_pderv[index2_k][index2_l]=cdl;
            
        }
    }
    
    return _SUCCESS_;
}

int lensing_addback_derv_cl_dl_tt(
                                  struct lensing * ple,
                                  double *cl_tt) {
    int index_l,index_k,l;
    
    for (index_l=0; index_l<ple->l_size; index_l++) {
        
        l = (int)ple->l_dl[index_l];
        
        for(index_k=0; index_k<ple->dl_size; index_k++){
            
            ple->cl_dl_tt_pderv[index_k][index_l] += cl_tt[l];
            
        }
    }
    return _SUCCESS_;
    
}

int lensing_delensed_derv_cl_te(
                                double **ksiX_dl_derv,
                                double *hla,
                                double *hla_P,
                                double **d20,
                                double *w8,
                                int nmu,
                                struct lensing * ple
                                ) {
    
    double clte_dl;
    int imu, index2_l,index2_k,l;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index2_l,index2_k,clte_dl) \
schedule (static)
    
    for(index2_k=0; index2_k<ple->dl_size; index2_k++){
        for(index2_l=0; index2_l<ple->dl_size; index2_l++){
            
            clte_dl=0.;
            
            for (imu=0;imu<nmu;imu++) {
                
                clte_dl += ksiX_dl_derv[index2_k][imu]*d20[imu][(int)ple->l_dl[index2_l]]*w8[imu]; /* loop could be optimized */
                
            }
            
            clte_dl *= 2.0*_PI_;
            
            ple->cl_dl_te_pderv[index2_k][index2_l]=clte_dl;
            
        }
    }
    
    return _SUCCESS_;
}

int lensing_addback_derv_cl_dl_te(
                                  struct lensing * ple,
                                  double *cl_te) {
    int index_l,index_k,l;
    
    for (index_l=0; index_l<ple->l_size; index_l++) {
        
        l = (int)ple->l_dl[index_l];
        
        for(index_k=0; index_k<ple->dl_size; index_k++){
            
            ple->cl_dl_te_pderv[index_k][index_l] += cl_te[l];
            
        }
    }
    return _SUCCESS_;
    
}

int lensing_delensed_derv_cl_ee_bb(
                                   double **ksip_dl_derv,
                                   double **ksim_dl_derv,
                                   double *hla,
                                   double **d22,
                                   double **d2m2,
                                   double *w8,
                                   int nmu,
                                   struct lensing * ple
                                   ) {
    
    double clp_dl, clm_dl, cl_ee, cl_bb;
    int imu;
    int index_l,index_k;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l,index_k,clp_dl,clm_dl,cl_ee,cl_bb) \
schedule (static)
    
    for(index_k=0; index_k < ple->dl_size; index_k++){
        
        for(index_l=0; index_l < ple->dl_size; index_l++){
            
            clp_dl=0; clm_dl=0;
            
            for (imu=0;imu<nmu;imu++) {
                
                clp_dl += ksip_dl_derv[index_k][imu]*d22[imu][(int)ple->l_dl[index_l]]*w8[imu];
                clm_dl += ksim_dl_derv[index_k][imu]*d2m2[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
            }
            
            cl_ee = (clp_dl+clm_dl)*_PI_;
            cl_bb = (clp_dl-clm_dl)*_PI_;
            
            ple->cl_dl_ee_pderv[index_k][index_l]=cl_ee;
            ple->cl_dl_bb_pderv[index_k][index_l]=cl_bb;
            
        }
        
    }
    
    return _SUCCESS_;
}

int lensing_addback_derv_cl_dl_ee_bb(
                                     struct lensing * ple,
                                     double * cl_ee,
                                     double * cl_bb) {
    
    int index_l, index_k, l;
    
    for (index_k=0; index_k<ple->dl_size; index_k++) {
        
        l = (int)ple->l_dl[index_l];
        
        for (index_l=0; index_l<ple->dl_size; index_l++) {
            
            ple->cl_dl_ee_pderv[index_k][index_l] += cl_ee[l];
            ple->cl_dl_bb_pderv[index_k][index_l] += cl_bb[l];
        }
    }
    return _SUCCESS_;
    
}


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
                           ) {
    
    int l,imu,index_l;
    double ll;
    ErrorMsg erreur;
    
    double N_EB;
    
    double *rat2, *rat3, *rat5, *rat6, *rat12, *rat13;
    double *fac31, *fac33, *fac11v2, *fac, *fac21v2, *fac32;
    
    class_alloc(rat2  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat3  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat5  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat6  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat12 ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat13 ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    
    class_alloc(fac31  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac33  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac11v2,(ple->l_lensed_max+1)*sizeof(double),erreur); // EB , EE, BB, TB, TE
    class_alloc(fac    ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac21v2,(ple->l_lensed_max+1)*sizeof(double),erreur); // TE
    class_alloc(fac32  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    
    //////// EE ////////
    double *zeta_e33,*zeta_e3m3,*zeta_e31,*zeta_e3m1;
    double *zeta_e11,*zeta_e1m1,*zeta_e22,*zeta_e2m2;
    double *zeta_e21,*zeta_e2m1,*zeta_e32,*zeta_e3m2;
    
    class_alloc(zeta_e33,nmu*sizeof(double),erreur);  class_alloc(zeta_e3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_e31,nmu*sizeof(double),erreur);  class_alloc(zeta_e3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_e11,nmu*sizeof(double),erreur);  class_alloc(zeta_e1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_e22,nmu*sizeof(double),erreur);  class_alloc(zeta_e2m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_e21,nmu*sizeof(double),erreur);  class_alloc(zeta_e2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_e32,nmu*sizeof(double),erreur);  class_alloc(zeta_e3m2,nmu*sizeof(double),erreur);
    
    //////// BB ////////
    double *zeta_b33,*zeta_b3m3,*zeta_b31,*zeta_b3m1;
    double *zeta_b11,*zeta_b1m1,*zeta_b22,*zeta_b2m2;
    double *zeta_b21,*zeta_b2m1,*zeta_b32,*zeta_b3m2;
    
    class_alloc(zeta_b33,nmu*sizeof(double),erreur);  class_alloc(zeta_b3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_b31,nmu*sizeof(double),erreur);  class_alloc(zeta_b3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_b11,nmu*sizeof(double),erreur);  class_alloc(zeta_b1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_b22,nmu*sizeof(double),erreur);  class_alloc(zeta_b2m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_b21,nmu*sizeof(double),erreur);  class_alloc(zeta_b2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_b32,nmu*sizeof(double),erreur);  class_alloc(zeta_b3m2,nmu*sizeof(double),erreur);
    
    
    for (l=2; l<=ple->l_lensed_max; l++) {
        
        ll = (double)l;
                
        rat2[l] = 1./(cl_ee_ln[l]+ple->pk_pn[l]);
        rat3[l] = 1./(cl_bb_ln[l]+ple->pk_pn[l]);
        rat5[l] = cl_ee_th[l]/(cl_ee_ln[l]+ple->pk_pn[l]);
        rat6[l] = cl_bb_th[l]/(cl_bb_ln[l]+ple->pk_pn[l]);
        rat12[l]= cl_ee_th[l]*cl_ee_th[l]/(cl_ee_ln[l]+ple->pk_pn[l]);
        rat13[l]= cl_bb_th[l]*cl_bb_th[l]/(cl_bb_ln[l]+ple->pk_pn[l]);
        
        fac[l] = 2.*ll+1.;
        fac11v2[l] = fac[l]*(ll+2.)*(ll-1.);
        fac21v2[l] = fac[l]*sqrt((ll+2.)*(ll-1.));
        fac31[l]   = fac[l]*sqrt((ll+2.)*(ll-1.)*(ll-2.)*(ll+3.));
        fac32[l]   = fac[l]*sqrt((ll-2.)*(ll+3.));
        fac33[l]   = fac[l]*(ll-2.)*(ll+3.);
        
        
        // Apply l_mask before computing reconstruction noise
        if (l<ple->recon_mask_lmin_E || l>ple->recon_mask_lmax_E) {
            rat2[l] = 0;
            rat5[l] = 0;
            rat12[l]= 0;
        }
        if (l<ple->recon_mask_lmin_B || l>ple->recon_mask_lmax_B) {
            rat3[l] = 0;
            rat6[l] = 0;
            rat13[l]= 0;
        }
    }


#pragma omp parallel for                        \
private (imu,l,ll) \
schedule (static)
    
    for (imu=0; imu<nmu; imu++) {
        
        //////// EE ////////
        zeta_e33[imu] = 0.; zeta_e3m3[imu] = 0.;
        zeta_e31[imu] = 0.; zeta_e3m1[imu] = 0.;
        zeta_e11[imu] = 0.; zeta_e1m1[imu] = 0.;
        zeta_e22[imu] = 0.; zeta_e2m2[imu] = 0.;
        zeta_e21[imu] = 0.; zeta_e2m1[imu] = 0.;
        zeta_e32[imu] = 0.; zeta_e3m2[imu] = 0.;
        
        //////// BB ////////
        zeta_b33[imu] = 0.; zeta_b3m3[imu] = 0.;
        zeta_b31[imu] = 0.; zeta_b3m1[imu] = 0.;
        zeta_b11[imu] = 0.; zeta_b1m1[imu] = 0.;
        zeta_b22[imu] = 0.; zeta_b2m2[imu] = 0.;
        zeta_b21[imu] = 0.; zeta_b2m1[imu] = 0.;
        zeta_b32[imu] = 0.; zeta_b3m2[imu] = 0.;
        
        for (l=2; l<=ple->l_lensed_max; l++) {
            
            ll = (double)l;
            
            //////// EE ////////
            zeta_e31[imu]  += rat12[l]*fac31[l]  *d31[imu][l];  zeta_e3m1[imu] += rat12[l]*fac31[l]  *d3m1[imu][l];
            zeta_e33[imu]  += rat12[l]*fac33[l]  *d33[imu][l];  zeta_e3m3[imu] += rat12[l]*fac33[l]  *d3m3[imu][l];
            zeta_e11[imu]  += rat12[l]*fac11v2[l]*d11[imu][l];  zeta_e1m1[imu] += rat12[l]*fac11v2[l]*d1m1[imu][l];
            zeta_e22[imu]  += rat2[l] *fac[l]    *d22[imu][l];  zeta_e2m2[imu] += rat2[l] *fac[l]    *d2m2[imu][l];
            zeta_e21[imu]  += rat5[l] *fac21v2[l]*d21[imu][l];  zeta_e2m1[imu] += rat5[l] *fac21v2[l]*d2m1[imu][l];
            zeta_e32[imu]  += rat5[l] *fac32[l]  *d32[imu][l];  zeta_e3m2[imu] += rat5[l] *fac32[l]  *d3m2[imu][l];
            
            //////// BB ////////
            zeta_b31[imu]  += rat13[l]*fac31[l]  *d31[imu][l];  zeta_b3m1[imu] += rat13[l]*fac31[l]  *d3m1[imu][l];
            zeta_b33[imu]  += rat13[l]*fac33[l]  *d33[imu][l];  zeta_b3m3[imu] += rat13[l]*fac33[l]  *d3m3[imu][l];
            zeta_b11[imu]  += rat13[l]*fac11v2[l]*d11[imu][l];  zeta_b1m1[imu] += rat13[l]*fac11v2[l]*d1m1[imu][l];
            zeta_b22[imu]  += rat3[l] *fac[l]    *d22[imu][l];  zeta_b2m2[imu] += rat3[l] *fac[l]    *d2m2[imu][l];
            zeta_b21[imu]  += rat6[l] *fac21v2[l]*d21[imu][l];  zeta_b2m1[imu] += rat6[l] *fac21v2[l]*d2m1[imu][l];
            zeta_b32[imu]  += rat6[l] *fac32[l]  *d32[imu][l];  zeta_b3m2[imu] += rat6[l] *fac32[l]  *d3m2[imu][l];
        }
        
        //////// EE ////////
        zeta_e33[imu]  /= 4.*_PI_; zeta_e3m3[imu]  /= 4.*_PI_;
        zeta_e31[imu]  /= 4.*_PI_; zeta_e3m1[imu]  /= 4.*_PI_;
        zeta_e11[imu]  /= 4.*_PI_; zeta_e1m1[imu]  /= 4.*_PI_;
        zeta_e22[imu]  /= 4.*_PI_; zeta_e2m2[imu]  /= 4.*_PI_;
        zeta_e21[imu]  /= 4.*_PI_; zeta_e2m1[imu]  /= 4.*_PI_;
        zeta_e32[imu]  /= 4.*_PI_; zeta_e3m2[imu]  /= 4.*_PI_;
        
        //////// BB ////////
        zeta_b33[imu]  /= 4.*_PI_; zeta_b3m3[imu]  /= 4.*_PI_;
        zeta_b31[imu]  /= 4.*_PI_; zeta_b3m1[imu]  /= 4.*_PI_;
        zeta_b11[imu]  /= 4.*_PI_; zeta_b1m1[imu]  /= 4.*_PI_;
        zeta_b22[imu]  /= 4.*_PI_; zeta_b2m2[imu]  /= 4.*_PI_;
        zeta_b21[imu]  /= 4.*_PI_; zeta_b2m1[imu]  /= 4.*_PI_;
        zeta_b32[imu]  /= 4.*_PI_; zeta_b3m2[imu]  /= 4.*_PI_;
    }

#pragma omp parallel for                        \
private (imu,index_l,ll,N_EB) \
schedule (static)
    
    for (index_l=0; index_l < ple->dl_size; index_l++)
    {
        ll = (double)ple->l_dl[index_l];
        
        N_EB = 0.;
        
        for (imu=0; imu<nmu; imu++) {
            
            ///////// N_EB;N_EB ///////////
            N_EB  +=  ( zeta_b22[imu]*zeta_e33[imu]
                       +zeta_b22[imu]*zeta_e11[imu]
                       -2.*zeta_b2m2[imu]*zeta_e3m1[imu]
                       +zeta_e22[imu]*zeta_b33[imu]
                       +zeta_e22[imu]*zeta_b11[imu]
                       -2.*zeta_e2m2[imu]*zeta_b3m1[imu]
                       +2.*( zeta_e3m2[imu]*zeta_b3m2[imu]
                            +zeta_e2m1[imu]*zeta_b2m1[imu]
                            -zeta_e32[imu]*zeta_b21[imu]
                            -zeta_b32[imu]*zeta_e21[imu])
                       )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
            
            N_EB -=   ( zeta_e2m2[imu]*zeta_b3m3[imu]
                       +zeta_e2m2[imu]*zeta_b1m1[imu]
                       -2.*zeta_e22[imu]*zeta_b31[imu]
                       +zeta_b2m2[imu]*zeta_e3m3[imu]
                       +zeta_b2m2[imu]*zeta_e1m1[imu]
                       -2.*zeta_b22[imu]*zeta_e31[imu]
                       -2.*( zeta_e32[imu]*zeta_b32[imu]
                            +zeta_e21[imu]*zeta_b21[imu]
                            +zeta_e3m2[imu]*zeta_b2m1[imu]
                            +zeta_b3m2[imu]*zeta_e2m1[imu])
                       )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
        }
        
        N_EB *= (ll*(ll+1.))*_PI_/4.;
        
        ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_eb] = 1.0/N_EB;
        
    }
    
    //////// EE ////////
    free(zeta_e33); free(zeta_e3m3);
    free(zeta_e31); free(zeta_e3m1);
    free(zeta_e11); free(zeta_e1m1);
    free(zeta_e22); free(zeta_e2m2);
    free(zeta_e21); free(zeta_e2m1);
    free(zeta_e32); free(zeta_e3m2);
    
    //////// BB ////////
    free(zeta_b33); free(zeta_b3m3);
    free(zeta_b31); free(zeta_b3m1);
    free(zeta_b11); free(zeta_b1m1);
    free(zeta_b22); free(zeta_b2m2);
    free(zeta_b21); free(zeta_b2m1);
    free(zeta_b32); free(zeta_b3m2);
    
    
    return _SUCCESS_;
    
}

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
                            double *cl_pp,   double *nl_ebeb,
                            struct lensing * ple,
                            struct harmonic * phr
                            ) {
    
    int l,imu,index_l;
    double ll;
    ErrorMsg erreur;
    
    double cl_bb_dl;
    
    double *rat1, *rat2, *rat3, *rat4;
    double *fac31, *fac33, *fac11v2, *fac, *fac11;
    
    class_alloc(rat1  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat2  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat3  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat4  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    
    class_alloc(fac31  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac33  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac11  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac11v2,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac    ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    
    //////// ratio 1 ////////
    double *zeta_a33,*zeta_a3m3,*zeta_a31,*zeta_a3m1;
    double *zeta_a11,*zeta_a1m1,*zeta_p1,*zeta_m1;
    
    class_alloc(zeta_a33,nmu*sizeof(double),erreur);  class_alloc(zeta_a3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_a31,nmu*sizeof(double),erreur);  class_alloc(zeta_a3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_a11,nmu*sizeof(double),erreur);  class_alloc(zeta_a1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_p1,nmu*sizeof(double),erreur);   class_alloc(zeta_m1,nmu*sizeof(double),erreur);
    
    //////// ratio 2 ////////
    double *zeta_b33,*zeta_b3m3,*zeta_b31,*zeta_b3m1;
    double *zeta_b11,*zeta_b1m1,*zeta_p2, *zeta_m2;
    
    class_alloc(zeta_b33,nmu*sizeof(double),erreur);  class_alloc(zeta_b3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_b31,nmu*sizeof(double),erreur);  class_alloc(zeta_b3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_b11,nmu*sizeof(double),erreur);  class_alloc(zeta_b1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_p2,nmu*sizeof(double),erreur);   class_alloc(zeta_m2,nmu*sizeof(double),erreur);
    
    for (l=2; l<=ple->l_lensed_max; l++) {
        
        ll = (double)l;
        
        rat1[l] = cl_ee_th[l];
        rat2[l] = cl_pp[l];
        
        rat3[l] = cl_ee_th[l]*cl_ee_ln[l]/(cl_ee_ln[l]+ple->pk_pn[l]);
        rat4[l] = cl_pp[l]*cl_pp[l]/(cl_pp[l]+nl_ebeb[l]);
        
        fac[l] = 2.*ll+1.;
        fac11[l]   = fac[l]*(ll+1.)*ll;
        fac11v2[l] = fac[l]*(ll+2.)*(ll-1.);
        fac31[l]   = fac[l]*sqrt((ll+2.)*(ll-1.)*(ll-2.)*(ll+3.));
        fac33[l]   = fac[l]*(ll-2.)*(ll+3.);
        
    }
    
    for (imu=0; imu<nmu; imu++) {
        
        //////// ratio 1 ////////
        zeta_a33[imu] = 0.; zeta_a3m3[imu] = 0.;
        zeta_a31[imu] = 0.; zeta_a3m1[imu] = 0.;
        zeta_a11[imu] = 0.; zeta_a1m1[imu] = 0.;
        zeta_p1[imu]= 0.;   zeta_m1[imu]= 0.;
        
        //////// ratio 2 ////////
        zeta_b33[imu] = 0.; zeta_b3m3[imu] = 0.;
        zeta_b31[imu] = 0.; zeta_b3m1[imu] = 0.;
        zeta_b11[imu] = 0.; zeta_b1m1[imu] = 0.;
        zeta_p2[imu]= 0.;   zeta_m2[imu]= 0.;
        
        for (l=2; l<=ple->l_lensed_max; l++) {
            
            ll = (double)l;
            
            //////// ratio 1 ////////
            zeta_a31[imu]  += rat1[l]*fac31[l]  *d31[imu][l];  zeta_a3m1[imu] += rat1[l]*fac31[l]  *d3m1[imu][l];
            zeta_a33[imu]  += rat1[l]*fac33[l]  *d33[imu][l];  zeta_a3m3[imu] += rat1[l]*fac33[l]  *d3m3[imu][l];
            zeta_a11[imu]  += rat1[l]*fac11v2[l]*d11[imu][l];  zeta_a1m1[imu] += rat1[l]*fac11v2[l]*d1m1[imu][l];
            zeta_p1[imu]   += rat2[l]*fac11[l]  *d11[imu][l];  zeta_m1[imu]   += rat2[l]*fac11[l]  *d1m1[imu][l];
            
            //////// ratio 2 ////////
            zeta_b31[imu]  += rat3[l]*fac31[l]  *d31[imu][l];  zeta_b3m1[imu] += rat3[l]*fac31[l]  *d3m1[imu][l];
            zeta_b33[imu]  += rat3[l]*fac33[l]  *d33[imu][l];  zeta_b3m3[imu] += rat3[l]*fac33[l]  *d3m3[imu][l];
            zeta_b11[imu]  += rat3[l]*fac11v2[l]*d11[imu][l];  zeta_b1m1[imu] += rat3[l]*fac11v2[l]*d1m1[imu][l];
            zeta_p2[imu]   += rat4[l]*fac11[l]  *d11[imu][l];  zeta_m2[imu]   += rat4[l]*fac11[l]  *d1m1[imu][l];
            
        }
        
        //////// ratio 1 ////////
        zeta_a33[imu]  /= 4.*_PI_; zeta_a3m3[imu]  /= 4.*_PI_;
        zeta_a31[imu]  /= 4.*_PI_; zeta_a3m1[imu]  /= 4.*_PI_;
        zeta_a11[imu]  /= 4.*_PI_; zeta_a1m1[imu]  /= 4.*_PI_;
        zeta_p1[imu]   /= 4.*_PI_; zeta_m1[imu]    /= 4.*_PI_;
        
        //////// ratio 2 ////////
        zeta_b33[imu]  /= 4.*_PI_; zeta_b3m3[imu]  /= 4.*_PI_;
        zeta_b31[imu]  /= 4.*_PI_; zeta_b3m1[imu]  /= 4.*_PI_;
        zeta_b11[imu]  /= 4.*_PI_; zeta_b1m1[imu]  /= 4.*_PI_;
        zeta_p2[imu]   /= 4.*_PI_; zeta_m2[imu]    /= 4.*_PI_;
        
    }
    
    for (index_l=0; index_l < ple->dl_size; index_l++)
    {
        ll = (double)ple->l_dl[index_l];
        
        cl_bb_dl = 0.;
        
        for (imu=0; imu<nmu; imu++) {
            
            ///////// cl_bb_lensed ///////////
            cl_bb_dl  +=  ( zeta_p1[imu]*zeta_a33[imu]
                           +zeta_p1[imu]*zeta_a11[imu]
                           +2.*zeta_m1[imu]*zeta_a3m1[imu]
                           )*d22[imu][(int)ple->l_dl[index_l]]*w8[imu];
            
            cl_bb_dl  -=  ( zeta_m1[imu]*zeta_a3m3[imu]
                           +zeta_m1[imu]*zeta_a1m1[imu]
                           +2.*zeta_p1[imu]*zeta_a31[imu]
                           )*d2m2[imu][(int)ple->l_dl[index_l]]*w8[imu];
        }
        
        cl_bb_dl *= _PI_/4.;
        
        ple->cl_delens[index_l*ple->dlt_size+ple->index_lt_dl_bb] = cl_bb_dl + ple->pk_pn[l];
        
    }
    
    //////// ratio 1 ////////
    free(zeta_a33); free(zeta_a3m3);
    free(zeta_a31); free(zeta_a3m1);
    free(zeta_a11); free(zeta_a1m1);
    free(zeta_p1);free(zeta_m1);
    
    //////// ratio 2 ////////
    free(zeta_b33); free(zeta_b3m3);
    free(zeta_b31); free(zeta_b3m1);
    free(zeta_b11); free(zeta_b1m1);
    free(zeta_p2);free(zeta_m2);
    
    
    return _SUCCESS_;
    
}


/**
 * This routine computes the minimum variance lensing noise reconstruction power spectra by Gaussian quadrature
 *
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param ple  Input/output: Pointer to the lensing structure
 * @return the error status
 */
int lensing_reconstr_nl_minvar(
                               double **d3m3,  double **d33,
                               double **d3m2,  double **d32,
                               double **d3m1,  double **d31,
                               double **d30,
                               double **d2m2,  double **d22,
                               double **d2m1,  double **d21,
                               double **d20,
                               double **d1m1,  double **d11,
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
                               ) {
    
    double *fac, *fac01, *fac11;
    double *fac01v2 = NULL, *fac11v2 = NULL, *fac11v3 = NULL;
    double *fac31   = NULL, *fac31v2 = NULL, *fac31v3 = NULL;
    double *fac33   = NULL, *fac32   = NULL, *fac20   = NULL;
    double *fac21   = NULL, *fac21v2 = NULL, *fac30   = NULL;
    
    int l,imu,index_l;
    double ll;
    
    ErrorMsg erreur;
    
    double *rat1, *rat2, *rat3, *rat4, *rat5, *rat6, *rat7, *rat8;
    double *rat9, *rat10,*rat11,*rat12,*rat13,*rat14,*rat15,*rat16;
    double *rat17,*rat18,*rat19,*rat20,*rat21,*rat22,*rat23,*rat24;
    double *rat25,*rat26,*rat27,*rat28,*rat29,*rat30,*rat31,*rat32;
    double *rat33;
    
    //////// TT ////////
    double *zeta_t00, *zeta_t01;
    double *zeta_t11,*zeta_t1m1;
    
    class_alloc(zeta_t00,nmu*sizeof(double),erreur);  class_alloc(zeta_t01,nmu*sizeof(double),erreur);
    class_alloc(zeta_t11,nmu*sizeof(double),erreur);  class_alloc(zeta_t1m1,nmu*sizeof(double),erreur);
    
    //////// TE ////////
    double *zeta_te31,*zeta_te3m1,*zeta_te11,*zeta_te1m1,*zeta_te21;
    double *zeta_te22,*zeta_te2m2,*zeta_te33,*zeta_te3m3,*zeta_te2m1;
    double *zeta_te01,*zeta_te30,*zeta_te00,*zeta_te11t,*zeta_te1m1t;
    
    class_alloc(zeta_te31,nmu*sizeof(double),erreur);  class_alloc(zeta_te3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_te11,nmu*sizeof(double),erreur);  class_alloc(zeta_te1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_te22,nmu*sizeof(double),erreur);  class_alloc(zeta_te2m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_te21,nmu*sizeof(double),erreur);  class_alloc(zeta_te2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_te33,nmu*sizeof(double),erreur);  class_alloc(zeta_te3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_te01,nmu*sizeof(double),erreur);  class_alloc(zeta_te30,nmu*sizeof(double),erreur);
    class_alloc(zeta_te11t,nmu*sizeof(double),erreur); class_alloc(zeta_te1m1t,nmu*sizeof(double),erreur);
    class_alloc(zeta_te00,nmu*sizeof(double),erreur);
    
    //////// EE ////////
    double *zeta_e33,*zeta_e3m3,*zeta_e31,*zeta_e3m1;
    double *zeta_e11,*zeta_e1m1,*zeta_e22,*zeta_e2m2;
    double *zeta_e21,*zeta_e2m1,*zeta_e32,*zeta_e3m2;
    
    
    class_alloc(zeta_e33,nmu*sizeof(double),erreur);  class_alloc(zeta_e3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_e31,nmu*sizeof(double),erreur);  class_alloc(zeta_e3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_e11,nmu*sizeof(double),erreur);  class_alloc(zeta_e1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_e22,nmu*sizeof(double),erreur);  class_alloc(zeta_e2m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_e21,nmu*sizeof(double),erreur);  class_alloc(zeta_e2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_e32,nmu*sizeof(double),erreur);  class_alloc(zeta_e3m2,nmu*sizeof(double),erreur);
    
    //////// BB ////////
    double *zeta_b33,*zeta_b3m3,*zeta_b31,*zeta_b3m1;
    double *zeta_b11,*zeta_b1m1,*zeta_b22,*zeta_b2m2;
    double *zeta_b21,*zeta_b2m1,*zeta_b32,*zeta_b3m2;
    
    class_alloc(zeta_b33,nmu*sizeof(double),erreur);  class_alloc(zeta_b3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_b31,nmu*sizeof(double),erreur);  class_alloc(zeta_b3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_b11,nmu*sizeof(double),erreur);  class_alloc(zeta_b1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_b22,nmu*sizeof(double),erreur);  class_alloc(zeta_b2m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_b21,nmu*sizeof(double),erreur);  class_alloc(zeta_b2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_b32,nmu*sizeof(double),erreur);  class_alloc(zeta_b3m2,nmu*sizeof(double),erreur);
    
    //////// TB ////////
    double *zeta_tb11,*zeta_tb1m1,*zeta_tb22,*zeta_tb2m2;
    double *zeta_tb31,*zeta_tb3m1,*zeta_tb33,*zeta_tb3m3;
    
    class_alloc(zeta_tb3m3,nmu*sizeof(double),erreur); class_alloc(zeta_tb11,nmu*sizeof(double),erreur);
    class_alloc(zeta_tb1m1,nmu*sizeof(double),erreur); class_alloc(zeta_tb22,nmu*sizeof(double),erreur);
    class_alloc(zeta_tb2m2,nmu*sizeof(double),erreur); class_alloc(zeta_tb31,nmu*sizeof(double),erreur);
    class_alloc(zeta_tb3m1,nmu*sizeof(double),erreur); class_alloc(zeta_tb33,nmu*sizeof(double),erreur);
    
    //////// TE;EE ////////
    double *zeta_1teee31, *zeta_1teee3m1, *zeta_1teee33, *zeta_1teee3m3, *zeta_1teee11, *zeta_1teee1m1;
    double *zeta_1teee22, *zeta_1teee2m2, *zeta_1teee32, *zeta_1teee3m2, *zeta_1teee21, *zeta_1teee2m1;
    double *zeta_2teee32, *zeta_2teee3m2, *zeta_2teee21, *zeta_2teee2m1, *zeta_1teee01, *zeta_1teee30;
    double *zeta_3teee21, *zeta_3teee2m1, *zeta_1teee20, *zeta_3teee31, *zeta_3teee3m1;
    double *zeta_3teee11, *zeta_3teee1m1;
    
    class_alloc(zeta_1teee31,nmu*sizeof(double),erreur);  class_alloc(zeta_1teee3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1teee22,nmu*sizeof(double),erreur);  class_alloc(zeta_1teee2m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_2teee32,nmu*sizeof(double),erreur);  class_alloc(zeta_2teee3m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_3teee21,nmu*sizeof(double),erreur);  class_alloc(zeta_3teee2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_3teee11,nmu*sizeof(double),erreur);  class_alloc(zeta_3teee1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1teee33,nmu*sizeof(double),erreur);  class_alloc(zeta_1teee3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_1teee32,nmu*sizeof(double),erreur);  class_alloc(zeta_1teee3m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_2teee21,nmu*sizeof(double),erreur);  class_alloc(zeta_2teee2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1teee20,nmu*sizeof(double),erreur);  class_alloc(zeta_3teee31,nmu*sizeof(double),erreur);
    class_alloc(zeta_1teee11,nmu*sizeof(double),erreur);  class_alloc(zeta_1teee1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1teee21,nmu*sizeof(double),erreur);  class_alloc(zeta_1teee2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1teee01,nmu*sizeof(double),erreur);  class_alloc(zeta_1teee30,nmu*sizeof(double),erreur);
    class_alloc(zeta_3teee3m1,nmu*sizeof(double),erreur);
    
    
    //////// TT;EE ////////
    double *zeta_1ttee20, *zeta_1ttee31, *zeta_1ttee11, *zeta_1ttee3m1, *zeta_1ttee1m1, *zeta_1ttee01;
    double *zeta_1ttee21, *zeta_1ttee30, *zeta_1ttee2m1;
    
    class_alloc(zeta_1ttee30,nmu*sizeof(double),erreur);  class_alloc(zeta_1ttee01,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttee21,nmu*sizeof(double),erreur);  class_alloc(zeta_1ttee2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttee31,nmu*sizeof(double),erreur);  class_alloc(zeta_1ttee3m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttee11,nmu*sizeof(double),erreur);  class_alloc(zeta_1ttee1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttee20,nmu*sizeof(double),erreur);
    
    //////// TT;TE ////////
    double *zeta_1ttte20, *zeta_1ttte31, *zeta_1ttte3m1, *zeta_1ttte11, *zeta_1ttte1m1;
    double *zeta_1ttte01, *zeta_2ttte01, *zeta_3ttte01,  *zeta_1ttte30, *zeta_1ttte21;
    double *zeta_1ttte2m1,*zeta_1ttte00, *zeta_2ttte11,  *zeta_2ttte1m1,*zeta_4ttte01;
    double *zeta_2ttte30, *zeta_2ttte21, *zeta_2ttte2m1, *zeta_2ttte00, *zeta_3ttte11;
    double *zeta_3ttte1m1,*zeta_2ttte20, *zeta_2ttte31,  *zeta_4ttte11, *zeta_2ttte3m1;
    double *zeta_4ttte1m1,*zeta_5ttte01, *zeta_6ttte01;
    
    class_alloc(zeta_1ttte20,nmu*sizeof(double),erreur);  class_alloc(zeta_2ttte11,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttte01,nmu*sizeof(double),erreur);  class_alloc(zeta_2ttte2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttte2m1,nmu*sizeof(double),erreur); class_alloc(zeta_3ttte11,nmu*sizeof(double),erreur);
    class_alloc(zeta_2ttte30,nmu*sizeof(double),erreur);  class_alloc(zeta_2ttte31,nmu*sizeof(double),erreur);
    class_alloc(zeta_3ttte1m1,nmu*sizeof(double),erreur); class_alloc(zeta_6ttte01,nmu*sizeof(double),erreur);
    class_alloc(zeta_4ttte1m1,nmu*sizeof(double),erreur); class_alloc(zeta_1ttte11,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttte31,nmu*sizeof(double),erreur);  class_alloc(zeta_1ttte30,nmu*sizeof(double),erreur);
    class_alloc(zeta_2ttte01,nmu*sizeof(double),erreur);  class_alloc(zeta_2ttte1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttte00,nmu*sizeof(double),erreur);  class_alloc(zeta_2ttte00,nmu*sizeof(double),erreur);
    class_alloc(zeta_2ttte21,nmu*sizeof(double),erreur);  class_alloc(zeta_4ttte11,nmu*sizeof(double),erreur);
    class_alloc(zeta_2ttte20,nmu*sizeof(double),erreur);  class_alloc(zeta_1ttte1m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_5ttte01,nmu*sizeof(double),erreur);  class_alloc(zeta_1ttte21,nmu*sizeof(double),erreur);
    class_alloc(zeta_1ttte3m1,nmu*sizeof(double),erreur); class_alloc(zeta_4ttte01,nmu*sizeof(double),erreur);
    class_alloc(zeta_3ttte01,nmu*sizeof(double),erreur);  class_alloc(zeta_2ttte3m1,nmu*sizeof(double),erreur);
    
    
    //////// TB;EB ////////
    double *zeta_1tbeb32, *zeta_1tbeb3m2, *zeta_1tbeb21, *zeta_1tbeb2m1, *zeta_2tbeb32, *zeta_2tbeb3m2;
    double *zeta_2tbeb21, *zeta_2tbeb2m1, *zeta_1tbeb33, *zeta_1tbeb3m3, *zeta_1tbeb11, *zeta_1tbeb1m1;
    double *zeta_1tbeb22, *zeta_1tbeb2m2, *zeta_1tbeb31, *zeta_1tbeb3m1;
    
    class_alloc(zeta_1tbeb32,nmu*sizeof(double),erreur);   class_alloc(zeta_1tbeb21,nmu*sizeof(double),erreur);
    class_alloc(zeta_2tbeb21,nmu*sizeof(double),erreur);   class_alloc(zeta_1tbeb33,nmu*sizeof(double),erreur);
    class_alloc(zeta_1tbeb22,nmu*sizeof(double),erreur);   class_alloc(zeta_1tbeb2m1,nmu*sizeof(double),erreur);
    class_alloc(zeta_1tbeb3m2,nmu*sizeof(double),erreur);  class_alloc(zeta_1tbeb3m3,nmu*sizeof(double),erreur);
    class_alloc(zeta_2tbeb2m1,nmu*sizeof(double),erreur);  class_alloc(zeta_2tbeb32,nmu*sizeof(double),erreur);
    class_alloc(zeta_1tbeb1m1,nmu*sizeof(double),erreur);  class_alloc(zeta_1tbeb11,nmu*sizeof(double),erreur);
    class_alloc(zeta_1tbeb2m2,nmu*sizeof(double),erreur);  class_alloc(zeta_2tbeb3m2,nmu*sizeof(double),erreur);
    class_alloc(zeta_1tbeb31,nmu*sizeof(double),erreur);   class_alloc(zeta_1tbeb3m1,nmu*sizeof(double),erreur);
    
    
    class_alloc(fac  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    
    class_alloc(fac01,(ple->l_lensed_max+1)*sizeof(double),erreur); // TT;TT - TT;TE
    class_alloc(fac11,(ple->l_lensed_max+1)*sizeof(double),erreur); // TT;TT - TT;TE
    
    class_alloc(fac01v2,(ple->l_lensed_max+1)*sizeof(double),erreur); // TE
    class_alloc(fac11v2,(ple->l_lensed_max+1)*sizeof(double),erreur); // EB , EE, BB, TB, TE
    class_alloc(fac11v3,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac31  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac31v2,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac31v3,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac33  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac32  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(fac20  ,(ple->l_lensed_max+1)*sizeof(double),erreur); // EE , EB , BB , TB
    class_alloc(fac21  ,(ple->l_lensed_max+1)*sizeof(double),erreur); // EE , EB , BB , TB
    class_alloc(fac21v2,(ple->l_lensed_max+1)*sizeof(double),erreur); // TE
    class_alloc(fac30  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    
    class_alloc(rat1  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat2  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat3  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat4  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat5  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat6  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat7  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat8  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat9  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat10  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat11  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat12  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat13  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat14  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat15  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat16  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat17  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat18  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat19  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat20  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat21  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat22  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat23  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat24  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat25  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat26  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat27  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat28  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat29  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat30  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat31  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat32  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    class_alloc(rat33  ,(ple->l_lensed_max+1)*sizeof(double),erreur);
    
    
    for (l=2; l<=ple->l_lensed_max; l++)
    {
        ll = (double)l;
        
        fac[l] = 2.*ll+1.;
        
        fac01[l]   = fac[l]*sqrt(ll*(ll+1.));
        fac01v2[l] = fac[l]*sqrt((ll+2.)*(ll-1.));
        
        fac11[l]   = fac[l]*ll*(ll+1.);
        fac11v2[l] = fac[l]*(ll+2.)*(ll-1.);
        fac11v3[l] = fac[l]*sqrt(ll*(ll+1.)*(ll+2.)*(ll-1.));
        
        fac20[l]   = fac[l]*sqrt(ll*(ll+1.)); // <- does this exist?
        fac21[l]   = fac[l]*sqrt(ll*(ll+1.));
        fac21v2[l] = fac[l]*sqrt((ll+2.)*(ll-1.));
        
        fac30[l]   = fac[l]*sqrt((ll-2.)*(ll+3.));
        
        fac31[l]   = fac[l]*sqrt((ll+2.)*(ll-1.)*(ll-2.)*(ll+3.));
        fac31v2[l] = fac[l]*sqrt(ll*(ll+1.)*(ll-2.)*(ll+3.));
        fac31v3[l] = fac[l]*sqrt((ll-2.)*(ll+3.)); // <- does this exist?
        
        fac32[l]   = fac[l]*sqrt((ll-2.)*(ll+3.));
        
        fac33[l]   = fac[l]*(ll-2.)*(ll+3.);
        
        
        rat1[l] = 1./(cl_tt_ln[l]+ple->pk_tn[l]);                                                                         // TT;TE:7
        rat2[l] = 1./(cl_ee_ln[l]+ple->pk_pn[l]);                                                                         //                 // TE;EE:1         //
        rat3[l] = 1./(cl_bb_ln[l]+ple->pk_pn[l]);                                                                         //
        rat4[l] = cl_tt_th[l]/(cl_tt_ln[l]+ple->pk_tn[l]);                                                                // TT;TE:4,
        rat33[l]= cl_tt_th[l]/(cl_ee_ln[l]+ple->pk_pn[l]);                                                                // TT;TE:10
        rat5[l] = cl_ee_th[l]/(cl_ee_ln[l]+ple->pk_pn[l]);                                                                //                 // TE;EE:9         //
        rat6[l] = cl_bb_th[l]/(cl_bb_ln[l]+ple->pk_pn[l]);                                                                //
        rat7[l] = cl_te_th[l]/(cl_tt_ln[l]+ple->pk_tn[l]);                                                                // TT;TE:5
        rat8[l] = cl_te_th[l]/(cl_ee_ln[l]+ple->pk_pn[l]);                                                                // TT;TE:15        // TE;EE:7,3       //
        rat9[l] = cl_tt_th[l]*cl_tt_th[l]/(cl_tt_ln[l]+ple->pk_tn[l]);                                                    //
        rat10[l]= cl_te_th[l]*cl_te_th[l]/(cl_tt_ln[l]+ple->pk_tn[l]);                                                    // TE
        rat11[l]= cl_te_th[l]*cl_te_th[l]/(cl_ee_ln[l]+ple->pk_pn[l]);                                                    //
        rat12[l]= cl_ee_th[l]*cl_ee_th[l]/(cl_ee_ln[l]+ple->pk_pn[l]);                                                    // EB
        rat13[l]= cl_bb_th[l]*cl_bb_th[l]/(cl_bb_ln[l]+ple->pk_pn[l]);                                                     //
        
        rat14[l] = cl_te_ln[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));                                 // TT;TE:1         // TE;EE:6,12      // TT;EE:13,12,7,2
        rat15[l] = cl_te_ln[l]*cl_te_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));                     // TT;TE:3         // TE;EE:10        //                 // TB;EB:3
        rat16[l] = cl_te_ln[l]*cl_tt_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));                     // TT;TE:6         //                 // TT;EE:16,9,6,3
        rat17[l] = cl_te_ln[l]*cl_te_th[l]*cl_tt_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));         // TT;TE:8
        rat18[l] = cl_te_ln[l]*cl_te_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_tt_ln[l]+ple->pk_tn[l]));                     // TT;TE:9
        rat19[l] = cl_te_ln[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_tt_ln[l]+ple->pk_tn[l]));                                 // TT;TE:11
        rat20[l] = 1./((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));                                          // TT;TE:13
        rat21[l] = cl_te_ln[l]*cl_tt_th[l]*cl_te_th[l]/(cl_tt_ln[l]+ple->pk_tn[l]);                                         // TT;TE:14
        rat22[l] = cl_te_th[l]*cl_tt_th[l]/(cl_ee_ln[l]+ple->pk_pn[l]);                                                    // TT;TE:12
        rat23[l] = cl_te_th[l]*cl_tt_th[l]/(cl_tt_ln[l]+ple->pk_tn[l]);                                                    // TT;TE:2
        rat24[l] = cl_te_ln[l]*cl_tt_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_tt_ln[l]+ple->pk_tn[l]));                     // TT;TE:16
        rat25[l] = cl_te_ln[l]*cl_te_th[l]*cl_ee_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));         //                 // TE;EE:2         // TT;EE:14,11,8,1
        rat26[l] = cl_te_th[l]*cl_ee_th[l]/((cl_ee_ln[l]+ple->pk_pn[l]));                                                 //                 // TE;EE:5,11      //
        rat27[l] = cl_te_ln[l]*cl_ee_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));                     //                 // TE;EE:8,4       // TT;EE:15,10,5,4
        rat28[l] = cl_te_ln[l]*cl_tt_th[l]*cl_ee_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));         //
        rat29[l] = cl_te_ln[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l])*(cl_tt_ln[l]+ple->pk_tn[l]));     //                 //                 //
        rat30[l] = (cl_bb_ln[l]+ple->pk_pn[l])/((cl_bb_ln[l]+ple->pk_pn[l])*(cl_bb_ln[l]+ple->pk_pn[l]));                 //                 //                 //                 // TB;EB:2
        rat31[l] = cl_te_ln[l]*cl_te_th[l]*cl_ee_th[l]/((cl_tt_ln[l]+ple->pk_tn[l])*(cl_ee_ln[l]+ple->pk_pn[l]));         //                 //                 //                 // TB;EB:1
        rat32[l] = (cl_bb_ln[l]+ple->pk_pn[l])*cl_bb_th[l]/((cl_bb_ln[l]+ple->pk_pn[l])*(cl_bb_ln[l]+ple->pk_pn[l]));     //                 //                 //                 // TB;EB:4
        

        // Apply l_mask before computing reconstruction noise
        if (l<ple->recon_mask_lmin_T || l>ple->recon_mask_lmax_T) {
            rat1[l] = 0;
            rat4[l] = 0;
            rat33[l]= 0;
            rat7[l] = 0;
            rat8[l] = 0;
            rat9[l] = 0;
            rat10[l]= 0;
            rat11[l]= 0;
            rat14[l] = 0;
            rat15[l] = 0;
            rat16[l] = 0;
            rat17[l] = 0;
            rat18[l] = 0;
            rat19[l] = 0;
            rat20[l] = 0;
            rat21[l] = 0;
            rat22[l] = 0;
            rat23[l] = 0;
            rat24[l] = 0;
            rat25[l] = 0;
            rat26[l] = 0;
            rat27[l] = 0;
            rat28[l] = 0;
            rat29[l] = 0;
            rat31[l] = 0;
        }
        if (l<ple->recon_mask_lmin_E || l>ple->recon_mask_lmax_E) {
            rat2[l] = 0;
            rat33[l]= 0;
            rat5[l] = 0;
            rat7[l] = 0;
            rat8[l] = 0;
            rat10[l]= 0;
            rat11[l]= 0;
            rat12[l]= 0;
            rat14[l] = 0;
            rat15[l] = 0;
            rat16[l] = 0;
            rat17[l] = 0;
            rat18[l] = 0;
            rat19[l] = 0;
            rat20[l] = 0;
            rat21[l] = 0;
            rat22[l] = 0;
            rat23[l] = 0;
            rat24[l] = 0;
            rat25[l] = 0;
            rat26[l] = 0;
            rat27[l] = 0;
            rat28[l] = 0;
            rat29[l] = 0;
            rat31[l] = 0;
        }
        if (l<ple->recon_mask_lmin_B || l>ple->recon_mask_lmax_B) {
            rat3[l] = 0;
            rat6[l] = 0;
            rat13[l]= 0;
            rat30[l]= 0;
            rat32[l]= 0;
        }
        
    }
    
    double N_TT;
    double N_TE;
    double N_EE;
    double N_BB;
    double N_EB;
    double N_TB;
    double N_TTTE;
    double N_TTEE;
    double N_TEEE;
    double N_TBEB;
    
    /** Integration by Gauss-Legendre quadrature. **/
#pragma omp parallel for                        \
private (imu,index_l,l,ll,N_TT,N_TE,N_EE,N_BB,N_EB,N_TB,    \
N_TTTE,N_TTEE,N_TEEE,N_TBEB) \
schedule (static)
    
    
    for (imu=0; imu<nmu; imu++) {
        
        //////// TT ////////
        zeta_t00[imu] = 0.; zeta_t01[imu]  = 0.;
        zeta_t11[imu] = 0.; zeta_t1m1[imu] = 0.;
        
        //////// TE ////////
        zeta_te31[imu] = 0.; zeta_te3m1[imu] = 0.;
        zeta_te11[imu] = 0.; zeta_te1m1[imu] = 0.;
        zeta_te22[imu] = 0.; zeta_te2m2[imu] = 0.;
        zeta_te21[imu] = 0.; zeta_te2m1[imu] = 0.;
        zeta_te33[imu] = 0.; zeta_te3m3[imu] = 0.;
        zeta_te01[imu] = 0.; zeta_te30[imu]  = 0.;
        zeta_te11t[imu]= 0.; zeta_te1m1t[imu]= 0.;
        zeta_te00[imu] = 0.;
        
        //////// EE ////////
        zeta_e33[imu] = 0.; zeta_e3m3[imu] = 0.;
        zeta_e31[imu] = 0.; zeta_e3m1[imu] = 0.;
        zeta_e11[imu] = 0.; zeta_e1m1[imu] = 0.;
        zeta_e22[imu] = 0.; zeta_e2m2[imu] = 0.;
        zeta_e21[imu] = 0.; zeta_e2m1[imu] = 0.;
        zeta_e32[imu] = 0.; zeta_e3m2[imu] = 0.;
        
        //////// BB ////////
        zeta_b33[imu] = 0.; zeta_b3m3[imu] = 0.;
        zeta_b31[imu] = 0.; zeta_b3m1[imu] = 0.;
        zeta_b11[imu] = 0.; zeta_b1m1[imu] = 0.;
        zeta_b22[imu] = 0.; zeta_b2m2[imu] = 0.;
        zeta_b21[imu] = 0.; zeta_b2m1[imu] = 0.;
        zeta_b32[imu] = 0.; zeta_b3m2[imu] = 0.;
        
        //////// TB ////////
        zeta_tb11[imu] = 0.; zeta_tb1m1[imu] = 0.;
        zeta_tb22[imu] = 0.; zeta_tb2m2[imu] = 0.;
        zeta_tb31[imu] = 0.; zeta_tb3m1[imu] = 0.;
        zeta_tb33[imu] = 0.; zeta_tb3m3[imu] = 0.;
        
        //////// TE;EE ////////
        zeta_1teee31[imu] = 0.;  zeta_1teee3m1[imu] = 0.;
        zeta_1teee22[imu] = 0.;  zeta_1teee2m2[imu] = 0.;
        zeta_2teee32[imu] = 0.;  zeta_2teee3m2[imu] = 0.;
        zeta_3teee21[imu] = 0.;  zeta_3teee2m1[imu] = 0.;
        zeta_3teee11[imu] = 0.;  zeta_3teee1m1[imu] = 0.;
        zeta_1teee33[imu] = 0.;  zeta_1teee3m3[imu] = 0.;
        zeta_1teee32[imu] = 0.;  zeta_1teee3m2[imu] = 0.;
        zeta_2teee21[imu] = 0.;  zeta_2teee2m1[imu] = 0.;
        zeta_1teee20[imu] = 0.;  zeta_3teee31[imu]  = 0.;
        zeta_1teee11[imu] = 0.;  zeta_1teee1m1[imu] = 0.;
        zeta_1teee21[imu] = 0.;  zeta_1teee2m1[imu] = 0.;
        zeta_1teee01[imu] = 0.;  zeta_1teee30[imu]  = 0.;
        zeta_3teee3m1[imu]= 0.;
        
        //////// TT;EE ////////
        zeta_1ttee21[imu] = 0.;  zeta_1ttee2m1[imu] = 0.;
        zeta_1ttee30[imu] = 0.;  zeta_1ttee01[imu]  = 0.;
        
        zeta_1ttee20[imu] = 0.;
        zeta_1ttee31[imu] = 0.;  zeta_1ttee3m1[imu] = 0.;
        zeta_1ttee11[imu] = 0.;  zeta_1ttee1m1[imu] = 0.;
        
        //////// TT;TE ////////
        zeta_1ttte20[imu]  = 0.;  zeta_2ttte11[imu]  = 0.;
        zeta_1ttte01[imu]  = 0.;  zeta_2ttte2m1[imu] = 0.;
        zeta_1ttte2m1[imu] = 0.;  zeta_3ttte11[imu]  = 0.;
        zeta_2ttte30[imu]  = 0.;  zeta_2ttte31[imu]  = 0.;
        zeta_3ttte1m1[imu] = 0.;  zeta_6ttte01[imu]  = 0.;
        zeta_4ttte1m1[imu] = 0.;  zeta_1ttte11[imu]  = 0.;
        zeta_1ttte31[imu]  = 0.;  zeta_1ttte30[imu]  = 0.;
        zeta_2ttte01[imu]  = 0.;  zeta_2ttte1m1[imu] = 0.;
        zeta_1ttte00[imu]  = 0.;  zeta_2ttte00[imu]  = 0.;
        zeta_2ttte21[imu]  = 0.;  zeta_4ttte11[imu]  = 0.;
        zeta_2ttte20[imu]  = 0.;  zeta_1ttte1m1[imu] = 0.;
        zeta_5ttte01[imu]  = 0.;  zeta_1ttte21[imu]  = 0.;
        zeta_1ttte3m1[imu] = 0.;  zeta_4ttte01[imu]  = 0.;
        zeta_3ttte01[imu]  = 0.;  zeta_2ttte3m1[imu] = 0.;
        
        //////// TB;EB ////////
        zeta_1tbeb32[imu]  = 0.;   zeta_1tbeb21[imu]  = 0.;
        zeta_2tbeb21[imu]  = 0.;   zeta_1tbeb33[imu]  = 0.;
        zeta_1tbeb22[imu]  = 0.;   zeta_1tbeb2m1[imu] = 0.;
        zeta_1tbeb3m2[imu] = 0.;  zeta_1tbeb3m3[imu]  = 0.;
        zeta_2tbeb2m1[imu] = 0.;  zeta_2tbeb32[imu]   = 0.;
        zeta_1tbeb1m1[imu] = 0.;  zeta_1tbeb11[imu]   = 0.;
        zeta_1tbeb2m2[imu] = 0.;  zeta_2tbeb3m2[imu]  = 0.;
        zeta_1tbeb31[imu]  = 0.;   zeta_1tbeb3m1[imu] = 0.;
        
        
        for (l=2; l<=ple->l_lensed_max; l++) {
            
            ll = (double)l;
            
            //////// TB ////////
            zeta_tb1m1[imu] += rat10[l]*fac11v2[l]*d1m1[imu][l];  zeta_tb11[imu] += rat10[l]*fac11v2[l]*d11[imu][l];
            zeta_tb3m1[imu] += rat10[l]*fac31[l]  *d3m1[imu][l];  zeta_tb31[imu] += rat10[l]*fac31[l]  *d31[imu][l];
            zeta_tb3m3[imu] += rat10[l]*fac33[l]  *d3m3[imu][l];  zeta_tb33[imu] += rat10[l]*fac33[l]  *d33[imu][l];
            zeta_tb2m2[imu] += rat3[l] *fac[l]    *d2m2[imu][l];  zeta_tb22[imu] += rat3[l] *fac[l]    *d22[imu][l];        //exists: zeta_b2m2
            
            //////// TT ////////
            zeta_t00[imu]  += rat1[l] *fac[l]   *d00[imu][l];                                                                //exists:zeta_te00
            zeta_t01[imu]  += rat4[l] *fac01[l] *d01[imu][l];
            zeta_t11[imu]  += rat9[l] *fac11[l] *d11[imu][l];
            zeta_t1m1[imu] += rat9[l] *fac11[l] *d1m1[imu][l];
            
            //////// TE ////////
            zeta_te31[imu] += rat10[l]*fac31[l]  *d31[imu][l];  zeta_te3m1[imu]  += rat10[l]*fac31[l]  *d3m1[imu][l];
            zeta_te11t[imu]+= rat10[l]*fac11v2[l]*d11[imu][l];  zeta_te1m1t[imu] += rat10[l]*fac11v2[l]*d1m1[imu][l];
            zeta_te33[imu] += rat10[l]*fac33[l]  *d33[imu][l];  zeta_te3m3[imu]  += rat10[l]*fac33[l]  *d3m3[imu][l];
            zeta_te11[imu] += rat11[l]*fac11v2[l]*d11[imu][l];  zeta_te1m1[imu]  += rat11[l]*fac11v2[l]*d1m1[imu][l];
            zeta_te22[imu] += rat2[l] *fac[l]    *d22[imu][l];  zeta_te2m2[imu]  += rat2[l] *fac[l]    *d2m2[imu][l];       //exists:zeta_e22
            zeta_te21[imu] += rat8[l] *fac21[l]  *d21[imu][l];  zeta_te2m1[imu]  += rat8[l] *fac21[l]  *d2m1[imu][l];
            zeta_te01[imu] += rat7[l] *fac01v2[l]*d01[imu][l];  zeta_te30[imu]   += rat7[l] *fac30[l]  *d30[imu][l];
            zeta_te00[imu] += fac[l]  *rat1[l]   *d00[imu][l];                                                                 //exists:zeta_t00
            
            //////// EE ////////
            zeta_e31[imu]  += rat12[l]*fac31[l]  *d31[imu][l];  zeta_e3m1[imu] += rat12[l]*fac31[l]  *d3m1[imu][l];
            zeta_e33[imu]  += rat12[l]*fac33[l]  *d33[imu][l];  zeta_e3m3[imu] += rat12[l]*fac33[l]  *d3m3[imu][l];
            zeta_e11[imu]  += rat12[l]*fac11v2[l]*d11[imu][l];  zeta_e1m1[imu] += rat12[l]*fac11v2[l]*d1m1[imu][l];
            zeta_e22[imu]  += rat2[l] *fac[l]    *d22[imu][l];  zeta_e2m2[imu] += rat2[l] *fac[l]    *d2m2[imu][l];         //exists:zeta_te22
            zeta_e21[imu]  += rat5[l] *fac21v2[l]*d21[imu][l];  zeta_e2m1[imu] += rat5[l] *fac21v2[l]*d2m1[imu][l];
            zeta_e32[imu]  += rat5[l] *fac32[l]  *d32[imu][l];  zeta_e3m2[imu] += rat5[l] *fac32[l]  *d3m2[imu][l];
            
            //////// BB ////////
            zeta_b31[imu]  += rat13[l]*fac31[l]  *d31[imu][l];  zeta_b3m1[imu] += rat13[l]*fac31[l]  *d3m1[imu][l];
            zeta_b33[imu]  += rat13[l]*fac33[l]  *d33[imu][l];  zeta_b3m3[imu] += rat13[l]*fac33[l]  *d3m3[imu][l];
            zeta_b11[imu]  += rat13[l]*fac11v2[l]*d11[imu][l];  zeta_b1m1[imu] += rat13[l]*fac11v2[l]*d1m1[imu][l];
            zeta_b22[imu]  += rat3[l] *fac[l]    *d22[imu][l];  zeta_b2m2[imu] += rat3[l] *fac[l]    *d2m2[imu][l];        //exists:zeta_tb2m2
            zeta_b21[imu]  += rat6[l] *fac21v2[l]*d21[imu][l];  zeta_b2m1[imu] += rat6[l] *fac21v2[l]*d2m1[imu][l];
            zeta_b32[imu]  += rat6[l] *fac32[l]  *d32[imu][l];  zeta_b3m2[imu] += rat6[l] *fac32[l]  *d3m2[imu][l];
            
            //////// TE;EE ////////
            zeta_1teee31[imu] += rat25[l] *fac31[l]  *d31[imu][l];  zeta_1teee3m1[imu] += rat25[l] *fac31[l]  *d3m1[imu][l]; //rat1->rat2
            zeta_1teee33[imu] += rat25[l] *fac33[l]  *d33[imu][l];  zeta_1teee3m3[imu] += rat25[l] *fac33[l]  *d3m3[imu][l];
            zeta_1teee11[imu] += rat25[l] *fac11v2[l]*d11[imu][l];  zeta_1teee1m1[imu] += rat25[l] *fac11v2[l]*d1m1[imu][l];
            zeta_1teee22[imu] += rat2[l]*fac[l]    *d22[imu][l];    zeta_1teee2m2[imu] += rat2[l]*fac[l]    *d2m2[imu][l]; //rat2->rat25
            
            zeta_1teee32[imu] += rat5[l] *fac32[l]  *d32[imu][l];  zeta_1teee3m2[imu] += rat5[l] *fac32[l]  *d3m2[imu][l]; // rat9
            zeta_1teee21[imu] += rat5[l] *fac21v2[l]*d21[imu][l];  zeta_1teee2m1[imu] += rat5[l] *fac21v2[l]*d2m1[imu][l]; // rat9
            
            zeta_2teee32[imu] += rat15[l]*fac32[l]  *d32[imu][l];  zeta_2teee3m2[imu] += rat15[l]*fac32[l]  *d3m2[imu][l]; // rat10
            zeta_2teee21[imu] += rat15[l]*fac21v2[l]*d21[imu][l];  zeta_2teee2m1[imu] += rat15[l]*fac21v2[l]*d2m1[imu][l]; // rat10
            
            zeta_1teee01[imu] += rat27[l]*fac01v2[l]*d01[imu][l];  zeta_1teee30[imu]  += rat27[l]*fac30[l]  *d30[imu][l];  // rat4
            zeta_3teee21[imu] += rat8[l] *fac21[l]  *d21[imu][l];  zeta_3teee2m1[imu] += rat8[l] *fac21[l]  *d2m1[imu][l]; // rat3
            
            zeta_1teee20[imu] += rat26[l]*fac[l]    *d20[imu][l];                                                           // rat12
            zeta_3teee31[imu] += rat14[l]*fac31v2[l]*d31[imu][l];  zeta_3teee3m1[imu] += rat14[l]*fac31v2[l]*d3m1[imu][l]; // rat11
            zeta_3teee11[imu] += rat14[l]*fac11v3[l]*d11[imu][l];  zeta_3teee1m1[imu] += rat14[l]*fac11v3[l]*d1m1[imu][l]; // rat11
            
            //////// TT;EE ////////
            zeta_1ttee01[imu] += rat27[l]*fac01v2[l]*d01[imu][l];  zeta_1ttee30[imu]  += rat27[l]*fac30[l]  *d30[imu][l];  //
            zeta_1ttee21[imu] += rat16[l]*fac21[l]  *d21[imu][l];  zeta_1ttee2m1[imu] += rat16[l]*fac21[l]  *d2m1[imu][l]; //
            
            zeta_1ttee20[imu] += rat14[l]*fac[l]    *d20[imu][l];                                                           //
            zeta_1ttee31[imu] += rat28[l]*fac31v2[l]*d31[imu][l];  zeta_1ttee3m1[imu] += rat28[l]*fac31v2[l]*d3m1[imu][l]; //
            zeta_1ttee11[imu] += rat28[l]*fac11v3[l]*d11[imu][l];  zeta_1ttee1m1[imu] += rat28[l]*fac11v3[l]*d1m1[imu][l]; //
            
            //////// TT;TE ////////
            zeta_1ttte20[imu] += rat14[l]*fac[l]    *d20[imu][l];                                                           // row1
            zeta_1ttte31[imu] += rat23[l]*fac31v2[l]*d31[imu][l];  zeta_1ttte3m1[imu] += rat23[l]*fac31v2[l]*d3m1[imu][l]; // row1
            zeta_1ttte11[imu] += rat23[l]*fac11v3[l]*d11[imu][l];  zeta_1ttte1m1[imu] += rat23[l]*fac11v3[l]*d1m1[imu][l]; // row1
            
            zeta_1ttte01[imu] += rat4[l] *fac01[l]  *d01[imu][l];  zeta_2ttte01[imu]  += rat15[l]*fac01[l]  *d01[imu][l];  // row1
            
            zeta_3ttte01[imu] += rat7[l] *fac01v2[l]*d01[imu][l];  zeta_1ttte30[imu]  += rat7[l] *fac30[l]  *d30[imu][l];  // row2
            zeta_1ttte21[imu] += rat16[l]*fac21[l]  *d21[imu][l];  zeta_1ttte2m1[imu] += rat16[l]*fac21[l]  *d2m1[imu][l]; // row2
            
            zeta_1ttte00[imu] += rat1[l] *fac[l]    *d00[imu][l];                                                          // row2  << probably same as some others
            zeta_2ttte11[imu] += rat17[l]*fac11[l]  *d11[imu][l];  zeta_2ttte1m1[imu] += rat17[l]*fac11[l]  *d1m1[imu][l]; // row2
            
            //////// TB;EB //////// // assumed some factors here - double check to make sure later.
            zeta_1tbeb32[imu] += rat32[l]*fac32[l]  *d32[imu][l];  zeta_1tbeb3m2[imu] += rat32[l]*fac32[l]  *d3m2[imu][l]; // row1
            zeta_1tbeb21[imu] += rat32[l]*fac21v2[l]*d21[imu][l];  zeta_1tbeb2m1[imu] += rat32[l]*fac21v2[l]*d2m1[imu][l]; // row1
            zeta_2tbeb32[imu] += rat15[l]*fac32[l]  *d32[imu][l];  zeta_2tbeb3m2[imu] += rat15[l]*fac32[l]  *d3m2[imu][l]; // row1
            zeta_2tbeb21[imu] += rat15[l]*fac21v2[l]*d21[imu][l];  zeta_2tbeb2m1[imu] += rat15[l]*fac21v2[l]*d2m1[imu][l]; // row1
            
            zeta_1tbeb33[imu] += rat31[l]*fac33[l]  *d33[imu][l];  zeta_1tbeb3m3[imu] += rat31[l]*fac33[l]  *d3m3[imu][l]; // row2
            zeta_1tbeb31[imu] += rat31[l]*fac31[l]  *d31[imu][l];  zeta_1tbeb3m1[imu] += rat31[l]*fac31[l]  *d3m1[imu][l]; // row2
            zeta_1tbeb11[imu] += rat31[l]*fac11v2[l]*d11[imu][l];  zeta_1tbeb1m1[imu] += rat31[l]*fac11v2[l]*d1m1[imu][l]; // row2
            zeta_1tbeb22[imu] += rat30[l]*fac[l]    *d22[imu][l];  zeta_1tbeb2m2[imu] += rat30[l]*fac[l]    *d2m2[imu][l]; // row2
            
        }
        
        //////// TT ////////
        zeta_t00[imu]  /= 4.*_PI_; zeta_t01[imu]   /= 4.*_PI_;
        zeta_t11[imu]  /= 4.*_PI_; zeta_t1m1[imu]  /= 4.*_PI_;
        
        //////// TE ////////
        zeta_te31[imu]  /= 4.*_PI_; zeta_te3m1[imu]  /= 4.*_PI_;
        zeta_te11[imu]  /= 4.*_PI_; zeta_te1m1[imu]  /= 4.*_PI_;
        zeta_te22[imu]  /= 4.*_PI_; zeta_te2m2[imu]  /= 4.*_PI_;
        zeta_te21[imu]  /= 4.*_PI_; zeta_te2m1[imu]  /= 4.*_PI_;
        zeta_te33[imu]  /= 4.*_PI_; zeta_te3m3[imu]  /= 4.*_PI_;
        zeta_te01[imu]  /= 4.*_PI_; zeta_te30[imu]   /= 4.*_PI_;
        zeta_te11t[imu] /= 4.*_PI_; zeta_te1m1t[imu] /= 4.*_PI_;
        zeta_te00[imu]  /= 4.*_PI_;
        
        //////// EE ////////
        zeta_e33[imu]  /= 4.*_PI_; zeta_e3m3[imu]  /= 4.*_PI_;
        zeta_e31[imu]  /= 4.*_PI_; zeta_e3m1[imu]  /= 4.*_PI_;
        zeta_e11[imu]  /= 4.*_PI_; zeta_e1m1[imu]  /= 4.*_PI_;
        zeta_e22[imu]  /= 4.*_PI_; zeta_e2m2[imu]  /= 4.*_PI_;
        zeta_e21[imu]  /= 4.*_PI_; zeta_e2m1[imu]  /= 4.*_PI_;
        zeta_e32[imu]  /= 4.*_PI_; zeta_e3m2[imu]  /= 4.*_PI_;
        
        //////// BB ////////
        zeta_b33[imu]  /= 4.*_PI_; zeta_b3m3[imu]  /= 4.*_PI_;
        zeta_b31[imu]  /= 4.*_PI_; zeta_b3m1[imu]  /= 4.*_PI_;
        zeta_b11[imu]  /= 4.*_PI_; zeta_b1m1[imu]  /= 4.*_PI_;
        zeta_b22[imu]  /= 4.*_PI_; zeta_b2m2[imu]  /= 4.*_PI_;
        zeta_b21[imu]  /= 4.*_PI_; zeta_b2m1[imu]  /= 4.*_PI_;
        zeta_b32[imu]  /= 4.*_PI_; zeta_b3m2[imu]  /= 4.*_PI_;
        
        //////// TB ////////
        zeta_tb11[imu]  /= 4.*_PI_; zeta_tb1m1[imu]  /= 4.*_PI_;
        zeta_tb22[imu]  /= 4.*_PI_; zeta_tb2m2[imu]  /= 4.*_PI_;
        zeta_tb31[imu]  /= 4.*_PI_; zeta_tb3m1[imu]  /= 4.*_PI_;
        zeta_tb33[imu]  /= 4.*_PI_; zeta_tb3m3[imu]  /= 4.*_PI_;
        
        //////// TE;EE ////////
        zeta_1teee22[imu]  /= 4.*_PI_;  zeta_1teee2m2[imu]  /= 4.*_PI_;
        zeta_1teee31[imu]  /= 4.*_PI_;  zeta_1teee3m1[imu]  /= 4.*_PI_;
        zeta_1teee11[imu]  /= 4.*_PI_;  zeta_1teee1m1[imu]  /= 4.*_PI_;
        zeta_1teee33[imu]  /= 4.*_PI_;  zeta_1teee3m3[imu]  /= 4.*_PI_;
        
        ////////
        zeta_1teee32[imu]  /= 4.*_PI_;  zeta_1teee3m2[imu]  /= 4.*_PI_;
        zeta_1teee21[imu]  /= 4.*_PI_;  zeta_1teee2m1[imu]  /= 4.*_PI_;
        
        zeta_2teee32[imu]  /= 4.*_PI_;  zeta_2teee3m2[imu]  /= 4.*_PI_;
        zeta_2teee21[imu]  /= 4.*_PI_;  zeta_2teee2m1[imu]  /= 4.*_PI_;
        ////////
        
        zeta_1teee01[imu]  /= 4.*_PI_;  zeta_1teee30[imu]   /= 4.*_PI_;
        zeta_3teee21[imu]  /= 4.*_PI_;  zeta_3teee2m1[imu]  /= 4.*_PI_;
        
        zeta_3teee31[imu]  /= 4.*_PI_;  zeta_3teee3m1[imu]  /= 4.*_PI_;
        zeta_3teee11[imu]  /= 4.*_PI_;  zeta_3teee1m1[imu]  /= 4.*_PI_;
        zeta_1teee20[imu]  /= 4.*_PI_;
        
        
        //////// TT;EE ////////
        
        zeta_1ttee21[imu]  /= 4.*_PI_;  zeta_1ttee2m1[imu]  /= 4.*_PI_;
        zeta_1ttee01[imu]  /= 4.*_PI_;  zeta_1ttee30[imu]   /= 4.*_PI_;
        
        zeta_1ttee31[imu]  /= 4.*_PI_;  zeta_1ttee3m1[imu]  /= 4.*_PI_;
        zeta_1ttee11[imu]  /= 4.*_PI_;  zeta_1ttee1m1[imu]  /= 4.*_PI_;
        zeta_1ttee20[imu]  /= 4.*_PI_;
        
        //////// TT;TE ////////
        
        ////////
        zeta_1ttte20[imu]  /= 4.*_PI_;
        zeta_1ttte31[imu]  /= 4.*_PI_;  zeta_1ttte3m1[imu] /= 4.*_PI_;
        zeta_1ttte11[imu]  /= 4.*_PI_;  zeta_1ttte1m1[imu] /= 4.*_PI_;
        
        zeta_1ttte01[imu]  /= 4.*_PI_;  zeta_2ttte01[imu]  /= 4.*_PI_;
        
        ////////
        zeta_1ttte30[imu]  /= 4.*_PI_;  zeta_3ttte01[imu]  /= 4.*_PI_;
        zeta_1ttte21[imu]  /= 4.*_PI_;  zeta_1ttte2m1[imu] /= 4.*_PI_;
        
        zeta_1ttte00[imu]  /= 4.*_PI_;
        zeta_2ttte11[imu]  /= 4.*_PI_;  zeta_2ttte1m1[imu] /= 4.*_PI_;
        
        ////////
        zeta_2ttte30[imu]  /= 4.*_PI_;  zeta_4ttte01[imu]  /= 4.*_PI_;
        zeta_2ttte21[imu]  /= 4.*_PI_;  zeta_2ttte2m1[imu] /= 4.*_PI_;
        
        zeta_2ttte00[imu]  /= 4.*_PI_;
        zeta_3ttte11[imu]  /= 4.*_PI_;  zeta_3ttte1m1[imu] /= 4.*_PI_;
        
        ////////
        zeta_2ttte20[imu]  /= 4.*_PI_;
        zeta_2ttte31[imu]  /= 4.*_PI_;  zeta_2ttte3m1[imu] /= 4.*_PI_;
        zeta_4ttte11[imu]  /= 4.*_PI_;  zeta_4ttte1m1[imu] /= 4.*_PI_;
        
        zeta_5ttte01[imu]  /= 4.*_PI_;  zeta_6ttte01[imu]  /= 4.*_PI_;
        
        
        //////// TB;EB ////////
        zeta_1tbeb32[imu]  /= 4.*_PI_;  zeta_1tbeb3m2[imu] /= 4.*_PI_;
        zeta_1tbeb21[imu]  /= 4.*_PI_;  zeta_1tbeb2m1[imu] /= 4.*_PI_;
        
        zeta_2tbeb21[imu]  /= 4.*_PI_;  zeta_2tbeb2m1[imu] /= 4.*_PI_;
        zeta_2tbeb32[imu]  /= 4.*_PI_;  zeta_2tbeb3m2[imu] /= 4.*_PI_;
        
        zeta_1tbeb22[imu]  /= 4.*_PI_;  zeta_1tbeb2m2[imu] /= 4.*_PI_;
        zeta_1tbeb33[imu]  /= 4.*_PI_;  zeta_1tbeb3m3[imu] /= 4.*_PI_;
        zeta_1tbeb11[imu]  /= 4.*_PI_;    zeta_1tbeb1m1[imu] /= 4.*_PI_;
        zeta_1tbeb31[imu]  /= 4.*_PI_;  zeta_1tbeb3m1[imu] /= 4.*_PI_;
    }


#pragma omp parallel for                        \
private (imu,index_l,l,ll,N_TT,N_TE,N_EE,N_BB,N_EB,N_TB,    \
N_TTTE,N_TTEE,N_TEEE,N_TBEB) \
schedule (static)
    
    for (index_l=0; index_l < ple->dl_size; index_l++)
    {
        ll = (double)ple->l_dl[index_l];
        
        N_TT = 0.;  N_TE = 0.;  N_EE = 0.;  N_BB = 0.;
        N_EB = 0.;  N_TB = 0.;  N_TTTE=0.;  N_TTEE=0.;
        N_TEEE=0.;  N_TBEB=0.;
        
        for (imu=0; imu<nmu; imu++) {
            
            ///////// N_TT;N_TT ///////////
            if(ple->has_nl_diag == _TRUE_ || ple->has_nl_eb_itr == _TRUE_)
            {
                
                N_TT += ( zeta_t00[imu]*zeta_t11[imu]
                         -zeta_t01[imu]*zeta_t01[imu]
                         )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                N_TT += ( zeta_t00[imu]*zeta_t1m1[imu]
                         +zeta_t01[imu]*zeta_t01[imu]
                         )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                ///////// N_EE;N_EE ///////////
                
                N_EE +=    ( zeta_e22[imu]*zeta_e33[imu]
                            +zeta_e22[imu]*zeta_e11[imu]
                            +2.*zeta_e2m2[imu]*zeta_e3m1[imu]
                            -( zeta_e3m2[imu]*zeta_e3m2[imu]
                              +zeta_e2m1[imu]*zeta_e2m1[imu]
                              +2.*zeta_e32[imu]*zeta_e21[imu])
                            )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                N_EE +=    ( zeta_e2m2[imu]*zeta_e3m3[imu]
                            +zeta_e2m2[imu]*zeta_e1m1[imu]
                            +2.*zeta_e22[imu]*zeta_e31[imu]
                            +( zeta_e32[imu]*zeta_e32[imu]
                              +zeta_e21[imu]*zeta_e21[imu]
                              -2.*zeta_e3m2[imu]*zeta_e2m1[imu])
                            )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                ///////// N_BB;N_BB ///////////
                
                N_BB +=    ( zeta_b22[imu]*zeta_b33[imu]
                            +zeta_b22[imu]*zeta_b11[imu]
                            +2.*zeta_b2m2[imu]*zeta_b3m1[imu]
                            -( zeta_b3m2[imu]*zeta_b3m2[imu]
                              +zeta_b2m1[imu]*zeta_b2m1[imu]
                              +2.*zeta_b32[imu]*zeta_b21[imu]) // typo!!!
                            )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                N_BB +=    ( zeta_b2m2[imu]*zeta_b3m3[imu]
                            +zeta_b2m2[imu]*zeta_b1m1[imu]
                            +2.*zeta_b22[imu]*zeta_b31[imu]
                            +( zeta_b32[imu]*zeta_b32[imu]
                              +zeta_b21[imu]*zeta_b21[imu]
                              -2.*zeta_b3m2[imu]*zeta_b2m1[imu])
                            )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                ///////// N_EB;N_EB ///////////
                if(ple->has_nl_eb_itr == _FALSE_)
                {
                    
                    N_EB  +=  ( zeta_b22[imu]*zeta_e33[imu]
                               +zeta_b22[imu]*zeta_e11[imu]
                               -2.*zeta_b2m2[imu]*zeta_e3m1[imu]
                               +zeta_e22[imu]*zeta_b33[imu]
                               +zeta_e22[imu]*zeta_b11[imu]
                               -2.*zeta_e2m2[imu]*zeta_b3m1[imu]
                               +2.*( zeta_e3m2[imu]*zeta_b3m2[imu]
                                    +zeta_e2m1[imu]*zeta_b2m1[imu]
                                    -zeta_e32[imu]*zeta_b21[imu]
                                    -zeta_b32[imu]*zeta_e21[imu])
                               )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                    
                    N_EB -=   ( zeta_e2m2[imu]*zeta_b3m3[imu]
                               +zeta_e2m2[imu]*zeta_b1m1[imu]
                               -2.*zeta_e22[imu]*zeta_b31[imu]
                               +zeta_b2m2[imu]*zeta_e3m3[imu]
                               +zeta_b2m2[imu]*zeta_e1m1[imu]
                               -2.*zeta_b22[imu]*zeta_e31[imu]
                               -2.*( zeta_e32[imu]*zeta_b32[imu]
                                    +zeta_e21[imu]*zeta_b21[imu]
                                    +zeta_e3m2[imu]*zeta_b2m1[imu]
                                    +zeta_b3m2[imu]*zeta_e2m1[imu])
                               )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                }
                ///////// N_TE;N_TE ///////////
                
                N_TE +=  (+( zeta_te22[imu]*zeta_te11t[imu]
                            +zeta_te22[imu]*zeta_te33[imu]
                            +2.*zeta_te2m2[imu]*zeta_te3m1[imu])
                          +4.*( zeta_te21[imu]*zeta_te01[imu]
                               -zeta_te2m1[imu]*zeta_te30[imu])
                          +4.*zeta_te00[imu]*zeta_te11[imu]
                          )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                N_TE +=  (+( zeta_te2m2[imu]*zeta_te1m1t[imu]
                            +zeta_te2m2[imu]*zeta_te3m3[imu]
                            +2.*zeta_te22[imu]*zeta_te31[imu])
                          +4.*( zeta_te2m1[imu]*zeta_te01[imu]
                               -zeta_te21[imu]*zeta_te30[imu])
                          +4.*zeta_te00[imu]*zeta_te1m1[imu]
                          )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                
                ///////// N_TB;N_TB ///////////
                
                N_TB += ( zeta_tb22[imu]*zeta_tb11[imu]
                         +zeta_tb22[imu]*zeta_tb33[imu]
                         -2.*zeta_tb2m2[imu]*zeta_tb3m1[imu]
                         )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                N_TB -= ( zeta_tb2m2[imu]*zeta_tb1m1[imu]
                         +zeta_tb2m2[imu]*zeta_tb3m3[imu]
                         -2.*zeta_tb22[imu]*zeta_tb31[imu]
                         )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                
                
                if(ple->has_nl_all == _TRUE_)
                {
                    
                    ///////// N_TT;N_TE ///////////
                    
                    N_TTTE += ( zeta_1ttte20[imu]*zeta_1ttte31[imu]    // row 1
                               +zeta_1ttte20[imu]*zeta_1ttte1m1[imu]
                               -2.*zeta_1ttte01[imu]*zeta_2ttte01[imu]
                               -zeta_1ttte2m1[imu]*zeta_1ttte30[imu]   // row 2
                               +zeta_1ttte21[imu] *zeta_3ttte01[imu]
                               +2.*zeta_1ttte00[imu]*zeta_2ttte11[imu]
                               )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                    
                    
                    N_TTTE += ( zeta_1ttte20[imu]*zeta_1ttte3m1[imu]    // row 1
                               +zeta_1ttte20[imu]*zeta_1ttte11[imu]
                               +2.*zeta_1ttte01[imu]*zeta_2ttte01[imu]
                               -zeta_1ttte21[imu]*zeta_1ttte30[imu]     // row 2
                               +zeta_1ttte2m1[imu] *zeta_3ttte01[imu]
                               +2.*zeta_1ttte00[imu]*zeta_2ttte1m1[imu]
                               )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                    
                    
                    ///////// N_TT;N_EE ///////////
                    
                    N_TTEE += ( zeta_1ttee20[imu]*zeta_1ttee31[imu]
                               +zeta_1ttee20[imu]*zeta_1ttee1m1[imu]
                               -zeta_1ttee2m1[imu]*zeta_1ttee30[imu]
                               +zeta_1ttee21[imu]*zeta_1ttee01[imu]
                               )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                    
                    
                    N_TTEE += ( zeta_1ttee20[imu]*zeta_1ttee3m1[imu]
                               +zeta_1ttee20[imu]*zeta_1ttee11[imu]
                               -zeta_1ttee21[imu]*zeta_1ttee30[imu]
                               +zeta_1ttee2m1[imu]*zeta_1ttee01[imu]
                               )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                    
                    
                    ///////// N_TE;N_EE ///////////
                    
                    N_TEEE +=  ( zeta_1teee22[imu]*zeta_1teee33[imu]     // row 1 + 5
                                +zeta_1teee22[imu]*zeta_1teee11[imu]
                                +2.*zeta_1teee2m2[imu]*zeta_1teee3m1[imu]
                                -zeta_1teee3m2[imu]*zeta_2teee3m2[imu]   // row 2 + 4
                                -zeta_1teee2m1[imu]*zeta_2teee2m1[imu]
                                -zeta_1teee32[imu]*zeta_2teee21[imu]
                                -zeta_2teee32[imu]*zeta_1teee21[imu]
                                +2.*(+zeta_1teee01[imu]*zeta_3teee21[imu]     // row 3 + 6
                                     -zeta_1teee30[imu]*zeta_3teee2m1[imu]
                                     +zeta_1teee20[imu]*zeta_3teee31[imu]
                                     +zeta_1teee20[imu]*zeta_3teee1m1[imu]
                                     )
                                )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                    
                    
                    N_TEEE +=  ( zeta_1teee2m2[imu]*zeta_1teee3m3[imu]   // row 1 + 5
                                +zeta_1teee2m2[imu]*zeta_1teee1m1[imu]
                                +2.*zeta_1teee22[imu]*zeta_1teee31[imu]
                                +zeta_1teee32[imu]*zeta_2teee32[imu]     // row 2 + 4
                                +zeta_1teee21[imu]*zeta_2teee21[imu]
                                -zeta_1teee3m2[imu]*zeta_2teee2m1[imu]
                                -zeta_2teee3m2[imu]*zeta_1teee2m1[imu]
                                +2.*(+zeta_1teee01[imu]*zeta_3teee2m1[imu]    // row 3 + 6
                                     -zeta_1teee30[imu]*zeta_3teee21[imu]
                                     +zeta_1teee20[imu]*zeta_3teee3m1[imu]
                                     +zeta_1teee20[imu]*zeta_3teee11[imu]
                                     )
                                )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                    
                    
                    
                    ///////// N_TB;N_EB ///////////
                    
                    N_TBEB += (-zeta_1tbeb3m2[imu]*zeta_2tbeb3m2[imu] // row 1
                               -zeta_1tbeb2m1[imu]*zeta_2tbeb2m1[imu]
                               +zeta_1tbeb32[imu]*zeta_2tbeb21[imu]
                               +zeta_2tbeb32[imu]*zeta_1tbeb21[imu]
                               +zeta_1tbeb22[imu]*zeta_1tbeb33[imu]   // row 2
                               +zeta_1tbeb22[imu]*zeta_1tbeb11[imu]
                               -2.*zeta_1tbeb2m2[imu]*zeta_1tbeb3m1[imu]
                               )*d11[imu][(int)ple->l_dl[index_l]]*w8[imu];
                    
                    N_TBEB -= (+zeta_1tbeb32[imu]*zeta_2tbeb32[imu]   // row 1
                               +zeta_1tbeb21[imu]*zeta_2tbeb21[imu]
                               +zeta_1tbeb3m2[imu]*zeta_2tbeb2m1[imu]
                               +zeta_2tbeb3m2[imu]*zeta_1tbeb2m1[imu]
                               +zeta_1tbeb2m2[imu]*zeta_1tbeb3m3[imu] // row 2
                               +zeta_1tbeb2m2[imu]*zeta_1tbeb1m1[imu]
                               -2.*zeta_1tbeb22[imu]*zeta_1tbeb31[imu]
                               )*d1m1[imu][(int)ple->l_dl[index_l]]*w8[imu];
                }
            }
        }
        
        if(ple->has_nl_diag == _TRUE_)
        {
            N_TT *= (ll*(ll+1.))*_PI_;
            
            ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_tt] = 1.0/N_TT;
            
            N_EE *= (ll*(ll+1.))*_PI_/4.;
            
            ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_ee] = 1.0/N_EE;
            
            N_BB *= (ll*(ll+1.))*_PI_/4.;
            
            ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_bb] = 1.0/N_BB;
            
            if(ple->has_nl_eb_itr == _FALSE_)
            {
                N_EB *= (ll*(ll+1.))*_PI_/4.;
                
                ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_eb] = 1.0/N_EB;
                
            }
            
            N_TE *= (ll*(ll+1.))*_PI_/4.;
            
            ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_te] = 1.0/N_TE;
            
            N_TB *= (ll*(ll+1.))*_PI_/4.;
            
            ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_tb] = 1.0/N_TB;
            
            if(ple->has_nl_all == _TRUE_)
            {
                N_TTTE *= (ll*(ll+1.))*_PI_/2.;
                
                ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_ttte] = N_TTTE/N_TT/N_TE;
                
                N_TTEE *= (ll*(ll+1.))*_PI_/2.;
                
                ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_ttee] = N_TTEE/N_EE/N_TT;
                
                N_TEEE *= (ll*(ll+1.))*_PI_/4.;
                
                ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_teee] = N_TEEE/N_TE/N_EE;
                
                N_TBEB *= (ll*(ll+1.))*_PI_/4.;
                
                ple->nl_rcn[index_l*ple->nlt_size+ple->index_nl_tbeb] = N_TBEB/N_TB/N_EB;
            }
        }
        
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //////// TT ////////
    free(zeta_t00); free(zeta_t01);
    free(zeta_t11); free(zeta_t1m1);
    
    //////// TE ////////
    free(zeta_te31); free(zeta_te3m1);
    free(zeta_te11); free(zeta_te1m1);
    free(zeta_te22); free(zeta_te2m2);
    free(zeta_te21); free(zeta_te2m1);
    free(zeta_te33); free(zeta_te3m3);
    free(zeta_te01); free(zeta_te30);
    free(zeta_te00);
    
    //////// EE ////////
    free(zeta_e33); free(zeta_e3m3);
    free(zeta_e31); free(zeta_e3m1);
    free(zeta_e11); free(zeta_e1m1);
    free(zeta_e22); free(zeta_e2m2);
    free(zeta_e21); free(zeta_e2m1);
    free(zeta_e32); free(zeta_e3m2);
    
    //////// BB ////////
    free(zeta_b33); free(zeta_b3m3);
    free(zeta_b31); free(zeta_b3m1);
    free(zeta_b11); free(zeta_b1m1);
    free(zeta_b22); free(zeta_b2m2);
    free(zeta_b21); free(zeta_b2m1);
    free(zeta_b32); free(zeta_b3m2);
    
    //////// TB ////////
    free(zeta_tb11); free(zeta_tb1m1);
    free(zeta_tb22); free(zeta_tb2m2);
    free(zeta_tb31); free(zeta_tb3m1);
    free(zeta_tb33); free(zeta_tb3m3);
    
    //////// TE;EE ////////
    free(zeta_1teee31);  free(zeta_1teee3m1);
    free(zeta_1teee22);  free(zeta_1teee2m2);
    free(zeta_2teee32);  free(zeta_2teee3m2);
    free(zeta_3teee21);  free(zeta_3teee2m1);
    free(zeta_3teee11);  free(zeta_3teee1m1);
    free(zeta_1teee33);  free(zeta_1teee3m3);
    free(zeta_1teee32);  free(zeta_1teee3m2);
    free(zeta_2teee21);  free(zeta_2teee2m1);
    free(zeta_1teee20);  free(zeta_3teee31);
    free(zeta_1teee11);  free(zeta_1teee1m1);
    free(zeta_1teee21);  free(zeta_1teee2m1);
    free(zeta_1teee01);  free(zeta_1teee30);
    free(zeta_3teee3m1);
    
    //////// TT;EE ////////
    free(zeta_1ttee20);
    free(zeta_1ttee31);  free(zeta_1ttee3m1);
    free(zeta_1ttee11);  free(zeta_1ttee1m1);
    
    free(zeta_1ttee30);  free(zeta_1ttee01);
    free(zeta_1ttee21);  free(zeta_1ttee2m1);
    
    //////// TT;TE ////////
    free(zeta_1ttte20);  free(zeta_2ttte11);
    free(zeta_1ttte01);  free(zeta_2ttte2m1);
    free(zeta_1ttte2m1); free(zeta_3ttte11);
    free(zeta_2ttte30);  free(zeta_2ttte31);
    free(zeta_3ttte1m1); free(zeta_6ttte01);
    free(zeta_4ttte1m1); free(zeta_1ttte11);
    free(zeta_1ttte31);  free(zeta_1ttte30);
    free(zeta_2ttte01);  free(zeta_2ttte1m1);
    free(zeta_1ttte00);  free(zeta_2ttte00);
    free(zeta_2ttte21);  free(zeta_4ttte11);
    free(zeta_2ttte20);  free(zeta_1ttte1m1);
    free(zeta_5ttte01);  free(zeta_1ttte21);
    free(zeta_1ttte3m1); free(zeta_4ttte01);
    free(zeta_3ttte01);  free(zeta_2ttte3m1);
    
    //////// TB;EB ////////
    free(zeta_1tbeb32);   free(zeta_1tbeb21);
    free(zeta_2tbeb21);   free(zeta_1tbeb33);
    free(zeta_1tbeb22);   free(zeta_1tbeb2m1);
    free(zeta_1tbeb3m2);  free(zeta_1tbeb3m3);
    free(zeta_2tbeb2m1);  free(zeta_2tbeb32);
    free(zeta_1tbeb1m1);  free(zeta_1tbeb11);
    free(zeta_1tbeb2m2);  free(zeta_2tbeb3m2);
    
    
    return _SUCCESS_;
}

/**
 * This routine computes the d00 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d00    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d00(
                double * mu,
                int num_mu,
                int lmax,
                double ** d00
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3;
    ErrorMsg erreur;
    
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    for (l=1; l<lmax; l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(2*ll+1)/(ll+1);
        fac2[l] = sqrt((2*ll+3)/(2*ll-1))*ll/(ll+1);
        fac3[l] = sqrt(2./(2*ll+3));
    }
    
#pragma omp parallel for                        \
private (index_mu,dlm1,dl,dlp1,l,ll)          \
schedule (static)
    
    for (index_mu=0;index_mu<num_mu;index_mu++) {
        dlm1=1.0/sqrt(2.); /* l=0 */
        d00[index_mu][0]=dlm1*sqrt(2.);
        dl=mu[index_mu] * sqrt(3./2.); /*l=1*/
        d00[index_mu][1]=dl*sqrt(2./3.);
        for(l=1;l<lmax;l++){
            ll=(double) l;
            /* sqrt((2l+1)/2)*d00 recurrence, supposed to be more stable */
            dlp1 = fac1[l]*mu[index_mu]*dl - fac2[l]*dlm1;
            d00[index_mu][l+1] = dlp1 * fac3[l];
            dlm1 = dl;
            dl = dlp1;
        }
    }
    free(fac1); free(fac2); free(fac3);
    return _SUCCESS_;
}


/*-----------------DLM--------------------*/
/*---------------------------------------*/
//
/**
 * This routine computes the d01 term <<--- SH: needed for lensing noise reconstruction.
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d01    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d01(
                double * mu,
                int num_mu,
                int lmax,
                double ** d01
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=2;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2.*ll+3.)*(2.*ll+1.)/(ll*(ll+2.)));
        fac3[l] = sqrt((2.*ll+3.)*(ll+1.)*(ll-1.)/(ll*(ll+2.)*(2.*ll-1)));
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
private (index_mu,dlm1,dl,dlp1,l,ll)          \
schedule (static)
    
    for (index_mu=0;index_mu<num_mu;index_mu++) {
        d01[index_mu][0]=0;
        dlm1=-sqrt(1.0-mu[index_mu])*sqrt(1.0+mu[index_mu])*sqrt(3./4.); /*l=1*/
        d01[index_mu][1]=dlm1 * sqrt(2./3.);
        dl=-(sqrt(1.0-mu[index_mu])*sqrt(1.0+mu[index_mu])*
             mu[index_mu]*sqrt(15./4.)); /*l=2*/
        d01[index_mu][2] = dl * sqrt(2./5.);
        for(l=2;l<lmax;l++){
            ll=(double) l;
            dlp1 = fac1[l]*mu[index_mu]*dl - fac3[l]*dlm1;
            d01[index_mu][l+1] = dlp1 * fac4[l];
            dlm1 = dl;
            dl = dlp1;
        }
    }
    free(fac1); free(fac3); free(fac4);
    return _SUCCESS_;
}

/**
 * This routine computes the d11 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d11    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d11(
                double * mu,
                int num_mu,
                int lmax,
                double ** d11
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=2;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/(ll*(ll+2));
        fac2[l] = 1.0/(ll*(ll+1.));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-1)*(ll+1)/(ll*(ll+2))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
private (index_mu,dlm1,dl,dlp1,l,ll)          \
schedule (static)
    
    for (index_mu=0;index_mu<num_mu;index_mu++) {
        d11[index_mu][0]=0;
        dlm1=(1.0+mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
        d11[index_mu][1]=dlm1 * sqrt(2./3.);
        dl=(1.0+mu[index_mu])/2.*(2.0*mu[index_mu]-1.0) * sqrt(5./2.); /*l=2*/
        d11[index_mu][2] = dl * sqrt(2./5.);
        for(l=2;l<lmax;l++){
            ll=(double) l;
            /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
            dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
            d11[index_mu][l+1] = dlp1 * fac4[l];
            dlm1 = dl;
            dl = dlp1;
        }
    }
    free(fac1); free(fac2); free(fac3); free(fac4);
    return _SUCCESS_;
}

/**
 * This routine computes the d1m1 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d1m1    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d1m1(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d1m1
                 ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=2;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/(ll*(ll+2));
        fac2[l] = 1.0/(ll*(ll+1.));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-1)*(ll+1)/(ll*(ll+2))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
private (index_mu,dlm1,dl,dlp1,l,ll)          \
schedule (static)
    
    for (index_mu=0;index_mu<num_mu;index_mu++) {
        d1m1[index_mu][0]=0;
        dlm1=(1.0-mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
        d1m1[index_mu][1]=dlm1 * sqrt(2./3.);
        dl=(1.0-mu[index_mu])/2.*(2.0*mu[index_mu]+1.0) * sqrt(5./2.); /*l=2*/
        d1m1[index_mu][2] = dl * sqrt(2./5.);
        for(l=2;l<lmax;l++){
            ll=(double) l;
            /* sqrt((2l+1)/2)*d1m1 recurrence, supposed to be more stable */
            dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
            d1m1[index_mu][l+1] = dlp1 * fac4[l];
            dlm1 = dl;
            dl = dlp1;
        }
    }
    free(fac1); free(fac2); free(fac3); free(fac4);
    return _SUCCESS_;
}

/*-----------------DLM--------------------*/
/*---------------------------------------*/
/**
 * This routine computes the d2m1 term <<--- SH: needed for lensing noise reconstruction.
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d2m1    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d2m1(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d2m1
                 ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=2;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-1)*(ll+3)*ll*(ll+2))) * (ll+1);
        fac2[l] = 2.0/(ll*(ll+1.));
        fac3[l] = sqrt( (2*ll+3)*(ll-2)*(ll+1)/((2*ll-1)*(ll+3)*ll))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
private (index_mu,dlm1,dl,dlp1,l,ll)          \
schedule (static)
    
    for (index_mu=0;index_mu<num_mu;index_mu++) {
        d2m1[index_mu][0]=0;
        dlm1=0; /*l=1*/
        d2m1[index_mu][1]=dlm1 * sqrt(2./3.);
        dl=0.5*sqrt(1.-mu[index_mu])*sqrt(1.+mu[index_mu])*(1.-mu[index_mu])*sqrt(5./2.); /*l=2*/
        d2m1[index_mu][2] = dl * sqrt(2./5.);
        for(l=2;l<lmax;l++){
            ll=(double) l;
            /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
            dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
            d2m1[index_mu][l+1] = dlp1 * fac4[l];
            dlm1 = dl;
            dl = dlp1;
        }
    }
    free(fac1); free(fac2); free(fac3); free(fac4);
    return _SUCCESS_;
}
/*---------------------------------------*/
/*---------------------------------------*/

/*-----------------DLM--------------------*/
/*---------------------------------------*/
/**
 * This routine computes the d21 term <<--- SH: needed for lensing noise reconstruction.
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d21    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d21(
                double * mu,
                int num_mu,
                int lmax,
                double ** d21
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=2;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-1)*(ll+3)*ll*(ll+2))) * (ll+1);
        fac2[l] = 2.0/(ll*(ll+1.));
        fac3[l] = sqrt( (2*ll+3)*(ll-2)*(ll+1)/((2*ll-1)*(ll+3)*ll))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
private (index_mu,dlm1,dl,dlp1,l,ll)          \
schedule (static)
    
    for (index_mu=0;index_mu<num_mu;index_mu++) {
        d21[index_mu][0]=0;
        dlm1=0; /*l=1*/
        d21[index_mu][1]=0;
        dl=0.5*sqrt((1.-mu[index_mu])*(1.+mu[index_mu]))*(1.+mu[index_mu])*sqrt(5./2.); /*l=2*/
        d21[index_mu][2] = dl * sqrt(2./5.);
        for(l=2;l<lmax;l++){
            ll=(double) l;
            /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
            dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
            d21[index_mu][l+1] = dlp1 * fac4[l];
            dlm1 = dl;
            dl = dlp1;
        }
    }
    free(fac1); free(fac2); free(fac3); free(fac4);
    return _SUCCESS_;
}
/*---------------------------------------*/
/*---------------------------------------*/



/**
 * This routine computes the d2m2 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d2m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d2m2(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d2m2
                 ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=2;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/((ll-1)*(ll+3));
        fac2[l] = 4.0/(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-2)*(ll+2)/((ll-1)*(ll+3))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
private (index_mu,dlm1,dl,dlp1,l,ll)          \
schedule (static)
    
    for (index_mu=0;index_mu<num_mu;index_mu++) {
        d2m2[index_mu][0]=0;
        dlm1=0.; /*l=1*/
        d2m2[index_mu][1]=0;
        dl=(1.0-mu[index_mu])*(1.0-mu[index_mu])/4. * sqrt(5./2.); /*l=2*/
        d2m2[index_mu][2] = dl * sqrt(2./5.);
        for(l=2;l<lmax;l++){
            ll=(double) l;
            /* sqrt((2l+1)/2)*d2m2 recurrence, supposed to be more stable */
            dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
            d2m2[index_mu][l+1] = dlp1 * fac4[l];
            dlm1 = dl;
            dl = dlp1;
        }
    }
    free(fac1); free(fac2); free(fac3); free(fac4);
    return _SUCCESS_;
}

/**
 * This routine computes the d22 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d22    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d22(
                double * mu,
                int num_mu,
                int lmax,
                double ** d22
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=2;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/((ll-1)*(ll+3));
        fac2[l] = 4.0/(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-2)*(ll+2)/((ll-1)*(ll+3))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
private (index_mu,dlm1,dl,dlp1,l,ll)          \
schedule (static)
    
    for (index_mu=0;index_mu<num_mu;index_mu++) {
        d22[index_mu][0]=0;
        dlm1=0.; /*l=1*/
        d22[index_mu][1]=0;
        dl=(1.0+mu[index_mu])*(1.0+mu[index_mu])/4. * sqrt(5./2.); /*l=2*/
        d22[index_mu][2] = dl * sqrt(2./5.);
        for(l=2;l<lmax;l++){
            ll=(double) l;
            /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
            dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
            d22[index_mu][l+1] = dlp1 * fac4[l];
            dlm1 = dl;
            dl = dlp1;
        }
    }
    free(fac1); free(fac2); free(fac3); free(fac4);
    return _SUCCESS_;
}

/**
 * This routine computes the d20 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d20    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d20(
                double * mu,
                int num_mu,
                int lmax,
                double ** d20
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=2;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-1)*(ll+3)));
        fac3[l] = sqrt((2*ll+3)*(ll-2)*(ll+2)/((2*ll-1)*(ll-1)*(ll+3)));
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    dlm1=1.0/sqrt(2.); /* l=0 */
    d00[index_mu][0]=dlm1*sqrt(2.);
    dl=mu[index_mu] * sqrt(3./2.); /*l=1*/
    d00[index_mu][1]=dl*sqrt(2./3.);
    for (l=1;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d00 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*mu[index_mu]*dl - fac2[l]*dlm1;
      d00[index_mu][l+1] = dlp1 * fac3[l];
      dlm1 = dl;
      dl = dlp1;
    }
    free(fac1); free(fac3); free(fac4);
    return _SUCCESS_;
}


/**
 ----------------DLM----------------------
 * This routine computes the d30 term  SH: <<<---------
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d30    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d30(
                double * mu,
                int num_mu,
                int lmax,
                double ** d30
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=3;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2.*ll+3.)*(2.*ll+1.)/((ll+4.)*(ll-2.)));
        fac3[l] = sqrt((2.*ll+3.)*(ll-3.)*(ll+3.)/((2.*ll-1.)*(ll+4.)*(ll-2.)));
        fac4[l] = sqrt(2./(2.*ll+3.));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d11[index_mu][0]=0;
    dlm1=(1.0+mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
    d11[index_mu][1]=dlm1 * sqrt(2./3.);
    dl=(1.0+mu[index_mu])/2.*(2.0*mu[index_mu]-1.0) * sqrt(5./2.); /*l=2*/
    d11[index_mu][2] = dl * sqrt(2./5.);
    for (l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
      d11[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d31 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d31    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d31(
                double * mu,
                int num_mu,
                int lmax,
                double ** d31
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=3;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-2)*(ll+4)*ll*(ll+2))) * (ll+1);
        fac2[l] = 3.0/(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1)*(ll-3)*(ll+3)*(ll-1)*(ll+1)/((ll-2)*(ll+4)*ll*(ll+2)))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d1m1[index_mu][0]=0;
    dlm1=(1.0-mu[index_mu])/2. * sqrt(3./2.); /*l=1*/
    d1m1[index_mu][1]=dlm1 * sqrt(2./3.);
    dl=(1.0-mu[index_mu])/2.*(2.0*mu[index_mu]+1.0) * sqrt(5./2.); /*l=2*/
    d1m1[index_mu][2] = dl * sqrt(2./5.);
    for (l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d1m1 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d1m1[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d3m1 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m1   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d3m1(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d3m1
                 ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=3;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-2)*(ll+4)*ll*(ll+2))) * (ll+1);
        fac2[l] = 3.0/(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1)*(ll-3)*(ll+3)*(ll-1)*(ll+1)/((ll-2)*(ll+4)*ll*(ll+2)))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d2m2[index_mu][0]=0;
    dlm1=0.; /*l=1*/
    d2m2[index_mu][1]=0;
    dl=(1.0-mu[index_mu])*(1.0-mu[index_mu])/4. * sqrt(5./2.); /*l=2*/
    d2m2[index_mu][2] = dl * sqrt(2./5.);
    for (l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d2m2 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d2m2[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d22 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d22    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d22(
                double * mu,
                int num_mu,
                int lmax,
                double ** d22
                ) {
  double ll, dlm1, dl, dlp1;
  int index_mu, l;
  double *fac1, *fac2, *fac3, *fac4;
  ErrorMsg erreur;
  class_alloc(fac1,lmax*sizeof(double),erreur);
  class_alloc(fac2,lmax*sizeof(double),erreur);
  class_alloc(fac3,lmax*sizeof(double),erreur);
  class_alloc(fac4,lmax*sizeof(double),erreur);
  for (l=2;l<lmax;l++) {
    ll = (double) l;
    fac1[l] = sqrt((2*ll+3)/(2*ll+1))*(ll+1)*(2*ll+1)/((ll-1)*(ll+3));
    fac2[l] = 4.0/(ll*(ll+1));
    fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-2)*(ll+2)/((ll-1)*(ll+3))*(ll+1)/ll;
    fac4[l] = sqrt(2./(2*ll+3));
  }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d22[index_mu][0]=0;
    dlm1=0.; /*l=1*/
    d22[index_mu][1]=0;
    dl=(1.0+mu[index_mu])*(1.0+mu[index_mu])/4. * sqrt(5./2.); /*l=2*/
    d22[index_mu][2] = dl * sqrt(2./5.);
    for (l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
      d22[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/*-----------------DLM--------------------*/
/*---------------------------------------*/
//
/**
 * This routine computes the d3m2 term <<--- SH: needed for lensing noise reconstruction.
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d3m2(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d3m2
                 ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=3;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-2)*(ll+4)*(ll-1)*(ll+3))) * (ll+1);
        fac2[l] = 6.0/(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1)*(ll-3)*(ll+2)/((ll+4)*(ll-1)))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d20[index_mu][0]=0;
    dlm1=0.; /*l=1*/
    d20[index_mu][1]=0;
    dl=sqrt(15.)/4.*(1-mu[index_mu]*mu[index_mu]); /*l=2*/
    d20[index_mu][2] = dl * sqrt(2./5.);
    for (l=2;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*mu[index_mu]*dl - fac3[l]*dlm1;
      d20[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac3); free(fac4);
  return _SUCCESS_;
}
/*---------------------------------------*/

/*-----------------DLM--------------------*/
/*---------------------------------------*/
//
/**
 * This routine computes the d32 term <<--- SH: needed for lensing noise reconstruction.
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d32   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d32(
                double * mu,
                int num_mu,
                int lmax,
                double ** d32
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=3;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-2)*(ll+4)*(ll-1)*(ll+3))) * (ll+1);
        fac2[l] = 6.0/(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1)*(ll-3)*(ll+2)/((ll+4)*(ll-1)))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d31[index_mu][0]=0;
    d31[index_mu][1]=0;
    dlm1=0.; /*l=2*/
    d31[index_mu][2]=0;
    dl=sqrt(105./2.)*(1+mu[index_mu])*(1+mu[index_mu])*(1-mu[index_mu])/8.; /*l=3*/
    d31[index_mu][3] = dl * sqrt(2./7.);
    for (l=3;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]-fac2[l])*dl - fac3[l]*dlm1;
      d31[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}
/*---------------------------------------*/


/**
 * This routine computes the d3m3 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m3   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d3m3(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d3m3
                 ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=3;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1))*(ll+1)/((ll-2)*(ll+4));
        fac2[l] = 9.0/(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-3)*(ll+3)*(l+1)/((ll-2)*(ll+4)*ll);
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d3m1[index_mu][0]=0;
    d3m1[index_mu][1]=0;
    dlm1=0.; /*l=2*/
    d3m1[index_mu][2]=0;
    dl=sqrt(105./2.)*(1+mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/8.; /*l=3*/
    d3m1[index_mu][3] = dl * sqrt(2./7.);
    for (l=3;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d3m1[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}


/*-----------------DLM-------------------*/
/*---------------------------------------*/
//
/**
 * This routine computes the d33 term <<--- SH: needed for lensing noise reconstruction.
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d33   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d33(
                double * mu,
                int num_mu,
                int lmax,
                double ** d33
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=3;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1))*(ll+1)/((ll-2)*(ll+4));
        fac2[l] = 9.0/(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-3)*(ll+3)*(l+1)/((ll-2)*(ll+4)*ll);
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d3m3[index_mu][0]=0;
    d3m3[index_mu][1]=0;
    dlm1=0.; /*l=2*/
    d3m3[index_mu][2]=0;
    dl=sqrt(7./2.)*(1-mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/8.; /*l=3*/
    d3m3[index_mu][3] = dl * sqrt(2./7.);
    for (l=3;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d3m3[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}
/*---------------------------------------*/
/*---------------------------------------*/




/**
 * This routine computes the d40 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d40    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d40(
                double * mu,
                int num_mu,
                int lmax,
                double ** d40
                ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=4;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-3)*(ll+5)));
        fac3[l] = sqrt((2*ll+3)*(ll-4)*(ll+4)/((2*ll-1)*(ll-3)*(ll+5)));
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d40[index_mu][0]=0;
    d40[index_mu][1]=0;
    d40[index_mu][2]=0;
    dlm1=0.; /*l=3*/
    d40[index_mu][3]=0;
    dl=sqrt(315.)*(1+mu[index_mu])*(1+mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/16.; /*l=4*/
    d40[index_mu][4] = dl * sqrt(2./9.);
    for (l=4;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*mu[index_mu]*dl - fac3[l]*dlm1;
      d40[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d4m2 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d4m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d4m2(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d4m2
                 ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=4;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1)/((ll-3)*(ll+5)*(ll-1)*(ll+3))) * (ll+1.);
        fac2[l] = 8./(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)*(ll-4)*(ll+4)*(ll-2)*(ll+2)/((2*ll-1)*(ll-3)*(ll+5)*(ll-1)*(ll+3)))*(ll+1)/ll;
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d4m2[index_mu][0]=0;
    d4m2[index_mu][1]=0;
    d4m2[index_mu][2]=0;
    dlm1=0.; /*l=3*/
    d4m2[index_mu][3]=0;
    dl=sqrt(126.)*(1+mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/16.; /*l=4*/
    d4m2[index_mu][4] = dl * sqrt(2./9.);
    for (l=4;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d4m2[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}

/**
 * This routine computes the d4m4 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d4m4   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int lensing_d4m4(
                 double * mu,
                 int num_mu,
                 int lmax,
                 double ** d4m4
                 ) {
    double ll, dlm1, dl, dlp1;
    int index_mu, l;
    double *fac1, *fac2, *fac3, *fac4;
    ErrorMsg erreur;
    class_alloc(fac1,lmax*sizeof(double),erreur);
    class_alloc(fac2,lmax*sizeof(double),erreur);
    class_alloc(fac3,lmax*sizeof(double),erreur);
    class_alloc(fac4,lmax*sizeof(double),erreur);
    for (l=4;l<lmax;l++) {
        ll = (double) l;
        fac1[l] = sqrt((2*ll+3)*(2*ll+1))*(ll+1)/((ll-3)*(ll+5));
        fac2[l] = 16./(ll*(ll+1));
        fac3[l] = sqrt((2*ll+3)/(2*ll-1))*(ll-4)*(ll+4)*(ll+1)/((ll-3)*(ll+5)*ll);
        fac4[l] = sqrt(2./(2*ll+3));
    }
#pragma omp parallel for                        \
  private (index_mu,dlm1,dl,dlp1,l,ll)          \
  schedule (static)

  for (index_mu=0;index_mu<num_mu;index_mu++) {
    d4m4[index_mu][0]=0;
    d4m4[index_mu][1]=0;
    d4m4[index_mu][2]=0;
    dlm1=0.; /*l=3*/
    d4m4[index_mu][3]=0;
    dl=sqrt(9./2.)*(1-mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])*(1-mu[index_mu])/16.; /*l=4*/
    d4m4[index_mu][4] = dl * sqrt(2./9.);
    for (l=4;l<lmax;l++){
      ll=(double) l;
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      dlp1 = fac1[l]*(mu[index_mu]+fac2[l])*dl - fac3[l]*dlm1;
      d4m4[index_mu][l+1] = dlp1 * fac4[l];
      dlm1 = dl;
      dl = dlp1;
    }
  }
  free(fac1); free(fac2); free(fac3); free(fac4);
  return _SUCCESS_;
}


