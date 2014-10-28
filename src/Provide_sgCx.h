#ifndef __Provide_sgCx__
#define __Provide_sgCx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgCx_val[EBINS][THBINS];
extern Double_t sgCx_err[EBINS][THBINS];
extern Double_t sgCx_unc[EBINS][THBINS];
extern Double_t sgCx_th[EBINS][THBINS];
extern Double_t sgCx_lo[EBINS];
extern Double_t sgCx_en[EBINS];
extern Double_t sgCx_hi[EBINS];
extern Double_t sgCx_wt[EBINS];
extern Double_t sgCx_sy[EBINS];
extern Double_t sgCx_sc[EBINS];
extern Bool_t sgCx_pre[EBINS];
extern Int_t sgCx_pts[EBINS];
extern Int_t sgCx_bin;
extern Char_t sgCx_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgCx(Char_t*, Double_t, Double_t);
void Sort_sgCx(Int_t l, Int_t r);
Int_t GetEnergyBins_sgCx(Int_t* bins=NULL);
Int_t GetNPts_sgCx();
Int_t ReadLine_sgCx(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgCx();
Double_t GetScale_sgCx();

//-----------------------------------------------------------------------------

#endif
