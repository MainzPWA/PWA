#ifndef __Provide_sgF__
#define __Provide_sgF__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgF_val[EBINS][THBINS];
extern Double_t sgF_err[EBINS][THBINS];
extern Double_t sgF_unc[EBINS][THBINS];
extern Double_t sgF_th[EBINS][THBINS];
extern Double_t sgF_lo[EBINS];
extern Double_t sgF_en[EBINS];
extern Double_t sgF_hi[EBINS];
extern Double_t sgF_wt[EBINS];
extern Double_t sgF_sy[EBINS];
extern Double_t sgF_sc[EBINS];
extern Bool_t sgF_pre[EBINS];
extern Int_t sgF_pts[EBINS];
extern Int_t sgF_bin;
extern Char_t sgF_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgF(Char_t*, Double_t, Double_t);
void Sort_sgF(Int_t l, Int_t r);
Int_t GetEnergyBins_sgF(Int_t* bins=NULL);
Int_t GetNPts_sgF();
Int_t ReadLine_sgF(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgF();
Double_t GetScale_sgF();

//-----------------------------------------------------------------------------

#endif
