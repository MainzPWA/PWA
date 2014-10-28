#ifndef __Provide_sgE__
#define __Provide_sgE__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgE_val[EBINS][THBINS];
extern Double_t sgE_err[EBINS][THBINS];
extern Double_t sgE_unc[EBINS][THBINS];
extern Double_t sgE_th[EBINS][THBINS];
extern Double_t sgE_lo[EBINS];
extern Double_t sgE_en[EBINS];
extern Double_t sgE_hi[EBINS];
extern Double_t sgE_wt[EBINS];
extern Double_t sgE_sy[EBINS];
extern Double_t sgE_sc[EBINS];
extern Bool_t sgE_pre[EBINS];
extern Int_t sgE_pts[EBINS];
extern Int_t sgE_bin;
extern Char_t sgE_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgE(Char_t*, Double_t, Double_t);
void Sort_sgE(Int_t l, Int_t r);
Int_t GetEnergyBins_sgE(Int_t* bins=NULL);
Int_t GetNPts_sgE();
Int_t ReadLine_sgE(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgE();
Double_t GetScale_sgE();

//-----------------------------------------------------------------------------

#endif
