#ifndef __Provide_sgS__
#define __Provide_sgS__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgS_val[EBINS][THBINS];
extern Double_t sgS_err[EBINS][THBINS];
extern Double_t sgS_unc[EBINS][THBINS];
extern Double_t sgS_th[EBINS][THBINS];
extern Double_t sgS_lo[EBINS];
extern Double_t sgS_en[EBINS];
extern Double_t sgS_hi[EBINS];
extern Double_t sgS_wt[EBINS];
extern Double_t sgS_sy[EBINS];
extern Double_t sgS_sc[EBINS];
extern Bool_t sgS_pre[EBINS];
extern Int_t sgS_pts[EBINS];
extern Int_t sgS_bin;
extern Char_t sgS_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgS(Char_t*, Double_t, Double_t);
void Sort_sgS(Int_t l, Int_t r);
Int_t GetEnergyBins_sgS(Int_t* bins=NULL);
Int_t GetNPts_sgS();
Int_t ReadLine_sgS(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgS();
Double_t GetScale_sgS();

//-----------------------------------------------------------------------------

#endif
