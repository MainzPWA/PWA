#ifndef __Provide_sgLx__
#define __Provide_sgLx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgLx_val[EBINS][THBINS];
extern Double_t sgLx_err[EBINS][THBINS];
extern Double_t sgLx_unc[EBINS][THBINS];
extern Double_t sgLx_th[EBINS][THBINS];
extern Double_t sgLx_lo[EBINS];
extern Double_t sgLx_en[EBINS];
extern Double_t sgLx_hi[EBINS];
extern Double_t sgLx_wt[EBINS];
extern Double_t sgLx_sy[EBINS];
extern Double_t sgLx_sc[EBINS];
extern Bool_t sgLx_pre[EBINS];
extern Int_t sgLx_pts[EBINS];
extern Int_t sgLx_bin;
extern Char_t sgLx_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgLx(Char_t*, Double_t, Double_t);
void Sort_sgLx(Int_t l, Int_t r);
Int_t GetEnergyBins_sgLx(Int_t* bins=NULL);
Int_t GetNPts_sgLx();
Int_t ReadLine_sgLx(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgLx();
Double_t GetScale_sgLx();

//-----------------------------------------------------------------------------

#endif
