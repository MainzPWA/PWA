#ifndef __Globals__
#define __Globals__

//-----------------------------------------------------------------------------

#include "Root.h"

//-----------------------------------------------------------------------------

extern Double_t gEnergy;

extern Double_t Fit_en[SOL][EBINS];
extern Double_t Fit_chi[SOL][EBINS];
extern Double_t Fit_pen[SOL][EBINS];
extern Double_t Fit_val[SOL][8*(LBINS-1)][EBINS];
extern Double_t Fit_err[SOL][8*(LBINS-1)][EBINS];
extern Int_t Fit_pts[SOL];
extern Double_t Model_en[EBINS];
extern Double_t Model_val[8*(LBINS-1)][EBINS];
extern Double_t Model_err[8*(LBINS-1)][EBINS];
extern Int_t Model_pts;

//Configuration options for fitting process
extern Bool_t FIX_EP[LBINS];
extern Bool_t FIX_EM[LBINS];
extern Bool_t FIX_MP[LBINS];
extern Bool_t FIX_MM[LBINS];
extern Bool_t FIX_RE_E0P;
extern Bool_t FIX_IM_E0P;
extern Bool_t ONLY_CROSS_S;
extern Bool_t ONLY_CROSS_F;
extern Bool_t SGT_ENERGIES;
extern Bool_t FIX_EP_PHASE[LBINS];
extern Bool_t FIX_EM_PHASE[LBINS];
extern Bool_t FIX_MP_PHASE[LBINS];
extern Bool_t FIX_MM_PHASE[LBINS];
extern Bool_t FIX_SCALES;
extern Bool_t PRINT_PENALTY;
extern Int_t ERROR_MODE;
extern Int_t PENALTY_MODE;
extern Int_t D_WAVES;
extern Int_t L_MAX;
extern Int_t ITERATIONS;
extern Int_t SOLUTIONS;
extern Double_t PENALTY[4];
extern Double_t WEIGHT[5];
extern Double_t SCALING;
extern Double_t VARIATION[2];
extern Double_t MIN_ENERGY;
extern Double_t MAX_ENERGY;
extern Double_t MASS_MESON;
extern Double_t MASS_INITIAL;
extern Double_t MASS_FINAL;
extern Double_t MASS2_MESON;
extern Double_t MASS2_INITIAL;
extern Double_t MASS2_FINAL;
extern Double_t THRESHOLD;
extern Double_t EPSILON;
extern Double_t EPSILON2;
extern Double_t BETA;

//-----------------------------------------------------------------------------

inline void Swap(Double_t* a, Double_t* b)
{
  Double_t tmp;

  tmp = *a;
  *a = *b;
  *b = tmp;
}

//-----------------------------------------------------------------------------

inline void Swap(Int_t* a, Int_t* b)
{
  Int_t tmp;

  tmp = *a;
  *a = *b;
  *b = tmp;
}

//-----------------------------------------------------------------------------

#endif
