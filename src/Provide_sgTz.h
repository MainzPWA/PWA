#ifndef __Provide_sgTz__
#define __Provide_sgTz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgTz_val[EBINS][THBINS];
extern Double_t sgTz_err[EBINS][THBINS];
extern Double_t sgTz_unc[EBINS][THBINS];
extern Double_t sgTz_th[EBINS][THBINS];
extern Double_t sgTz_lo[EBINS];
extern Double_t sgTz_en[EBINS];
extern Double_t sgTz_hi[EBINS];
extern Double_t sgTz_wt[EBINS];
extern Double_t sgTz_sy[EBINS];
extern Double_t sgTz_sc[EBINS];
extern Bool_t sgTz_pre[EBINS];
extern Int_t sgTz_pts[EBINS];
extern Int_t sgTz_bin;
extern Char_t sgTz_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgTz(Char_t*, Double_t, Double_t);
void Sort_sgTz(Int_t l, Int_t r);
Int_t GetEnergyBins_sgTz(Int_t* bins=NULL);
Int_t GetNPts_sgTz();
Int_t ReadLine_sgTz(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgTz();
Double_t GetScale_sgTz();

//-----------------------------------------------------------------------------

#endif
