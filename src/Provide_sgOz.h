#ifndef __Provide_sgOz__
#define __Provide_sgOz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgOz_val[EBINS][THBINS];
extern Double_t sgOz_err[EBINS][THBINS];
extern Double_t sgOz_unc[EBINS][THBINS];
extern Double_t sgOz_th[EBINS][THBINS];
extern Double_t sgOz_lo[EBINS];
extern Double_t sgOz_en[EBINS];
extern Double_t sgOz_hi[EBINS];
extern Double_t sgOz_wt[EBINS];
extern Double_t sgOz_sy[EBINS];
extern Double_t sgOz_sc[EBINS];
extern Bool_t sgOz_pre[EBINS];
extern Int_t sgOz_pts[EBINS];
extern Int_t sgOz_bin;
extern Char_t sgOz_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgOz(Char_t*, Double_t, Double_t);
void Sort_sgOz(Int_t l, Int_t r);
Int_t GetEnergyBins_sgOz(Int_t* bins=NULL);
Int_t GetNPts_sgOz();
Int_t ReadLine_sgOz(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgOz();
Double_t GetScale_sgOz();

//-----------------------------------------------------------------------------

#endif
