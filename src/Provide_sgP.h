#ifndef __Provide_sgP__
#define __Provide_sgP__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgP_val[EBINS][THBINS];
extern Double_t sgP_err[EBINS][THBINS];
extern Double_t sgP_unc[EBINS][THBINS];
extern Double_t sgP_th[EBINS][THBINS];
extern Double_t sgP_lo[EBINS];
extern Double_t sgP_en[EBINS];
extern Double_t sgP_hi[EBINS];
extern Double_t sgP_wt[EBINS];
extern Double_t sgP_sy[EBINS];
extern Double_t sgP_sc[EBINS];
extern Bool_t sgP_pre[EBINS];
extern Int_t sgP_pts[EBINS];
extern Int_t sgP_bin;
extern Char_t sgP_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgP(Char_t*, Double_t, Double_t);
void Sort_sgP(Int_t l, Int_t r);
Int_t GetEnergyBins_sgP(Int_t* bins=NULL);
Int_t GetNPts_sgP();
Int_t ReadLine_sgP(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgP();
Double_t GetScale_sgP();

//-----------------------------------------------------------------------------

#endif
