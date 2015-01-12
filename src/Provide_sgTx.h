#ifndef __Provide_sgTx__
#define __Provide_sgTx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgTx_val[EBINS][THBINS];
extern Double_t sgTx_err[EBINS][THBINS];
extern Double_t sgTx_unc[EBINS][THBINS];
extern Double_t sgTx_th[EBINS][THBINS];
extern Double_t sgTx_lo[EBINS];
extern Double_t sgTx_en[EBINS];
extern Double_t sgTx_hi[EBINS];
extern Double_t sgTx_wt[EBINS];
extern Double_t sgTx_sy[EBINS];
extern Double_t sgTx_sc[EBINS];
extern Bool_t sgTx_pre[EBINS];
extern Int_t sgTx_pts[EBINS];
extern Int_t sgTx_bin;
extern Char_t sgTx_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgTx(Char_t*, Double_t, Double_t);
void Sort_sgTx(Int_t l, Int_t r);
Int_t GetEnergyBins_sgTx(Int_t* bins=NULL);
Int_t GetNPts_sgTx();
Int_t ReadLine_sgTx(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgTx();
Double_t GetScale_sgTx();

//-----------------------------------------------------------------------------

#endif
