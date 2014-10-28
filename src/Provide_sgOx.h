#ifndef __Provide_sgOx__
#define __Provide_sgOx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgOx_val[EBINS][THBINS];
extern Double_t sgOx_err[EBINS][THBINS];
extern Double_t sgOx_unc[EBINS][THBINS];
extern Double_t sgOx_th[EBINS][THBINS];
extern Double_t sgOx_lo[EBINS];
extern Double_t sgOx_en[EBINS];
extern Double_t sgOx_hi[EBINS];
extern Double_t sgOx_wt[EBINS];
extern Double_t sgOx_sy[EBINS];
extern Double_t sgOx_sc[EBINS];
extern Bool_t sgOx_pre[EBINS];
extern Int_t sgOx_pts[EBINS];
extern Int_t sgOx_bin;
extern Char_t sgOx_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgOx(Char_t*, Double_t, Double_t);
void Sort_sgOx(Int_t l, Int_t r);
Int_t GetEnergyBins_sgOx(Int_t* bins=NULL);
Int_t GetNPts_sgOx();
Int_t ReadLine_sgOx(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgOx();
Double_t GetScale_sgOx();

//-----------------------------------------------------------------------------

#endif
