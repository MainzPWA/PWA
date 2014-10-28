#ifndef __Provide_sgCz__
#define __Provide_sgCz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgCz_val[EBINS][THBINS];
extern Double_t sgCz_err[EBINS][THBINS];
extern Double_t sgCz_unc[EBINS][THBINS];
extern Double_t sgCz_th[EBINS][THBINS];
extern Double_t sgCz_lo[EBINS];
extern Double_t sgCz_en[EBINS];
extern Double_t sgCz_hi[EBINS];
extern Double_t sgCz_wt[EBINS];
extern Double_t sgCz_sy[EBINS];
extern Double_t sgCz_sc[EBINS];
extern Bool_t sgCz_pre[EBINS];
extern Int_t sgCz_pts[EBINS];
extern Int_t sgCz_bin;
extern Char_t sgCz_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgCz(Char_t*, Double_t, Double_t);
void Sort_sgCz(Int_t l, Int_t r);
Int_t GetEnergyBins_sgCz(Int_t* bins=NULL);
Int_t GetNPts_sgCz();
Int_t ReadLine_sgCz(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgCz();
Double_t GetScale_sgCz();

//-----------------------------------------------------------------------------

#endif
