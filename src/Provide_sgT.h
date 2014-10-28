#ifndef __Provide_sgT__
#define __Provide_sgT__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgT_val[EBINS][THBINS];
extern Double_t sgT_err[EBINS][THBINS];
extern Double_t sgT_unc[EBINS][THBINS];
extern Double_t sgT_th[EBINS][THBINS];
extern Double_t sgT_lo[EBINS];
extern Double_t sgT_en[EBINS];
extern Double_t sgT_hi[EBINS];
extern Double_t sgT_wt[EBINS];
extern Double_t sgT_sy[EBINS];
extern Double_t sgT_sc[EBINS];
extern Bool_t sgT_pre[EBINS];
extern Int_t sgT_pts[EBINS];
extern Int_t sgT_bin;
extern Char_t sgT_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgT(Char_t*, Double_t, Double_t);
void Sort_sgT(Int_t l, Int_t r);
Int_t GetEnergyBins_sgT(Int_t* bins=NULL);
Int_t GetNPts_sgT();
Int_t ReadLine_sgT(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgT();
Double_t GetScale_sgT();

//-----------------------------------------------------------------------------

#endif
