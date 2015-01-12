#ifndef __Provide_sgLz__
#define __Provide_sgLz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgLz_val[EBINS][THBINS];
extern Double_t sgLz_err[EBINS][THBINS];
extern Double_t sgLz_unc[EBINS][THBINS];
extern Double_t sgLz_th[EBINS][THBINS];
extern Double_t sgLz_lo[EBINS];
extern Double_t sgLz_en[EBINS];
extern Double_t sgLz_hi[EBINS];
extern Double_t sgLz_wt[EBINS];
extern Double_t sgLz_sy[EBINS];
extern Double_t sgLz_sc[EBINS];
extern Bool_t sgLz_pre[EBINS];
extern Int_t sgLz_pts[EBINS];
extern Int_t sgLz_bin;
extern Char_t sgLz_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgLz(Char_t*, Double_t, Double_t);
void Sort_sgLz(Int_t l, Int_t r);
Int_t GetEnergyBins_sgLz(Int_t* bins=NULL);
Int_t GetNPts_sgLz();
Int_t ReadLine_sgLz(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgLz();
Double_t GetScale_sgLz();

//-----------------------------------------------------------------------------

#endif
