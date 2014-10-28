#ifndef __Provide_S__
#define __Provide_S__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t S_val[EBINS][THBINS];
extern Double_t S_err[EBINS][THBINS];
extern Double_t S_unc[EBINS][THBINS];
extern Double_t S_th[EBINS][THBINS];
extern Double_t S_lo[EBINS];
extern Double_t S_en[EBINS];
extern Double_t S_hi[EBINS];
extern Double_t S_wt[EBINS];
extern Double_t S_sy[EBINS];
extern Double_t S_sc[EBINS];
extern Bool_t S_pre[EBINS];
extern Int_t S_pts[EBINS];
extern Int_t S_bin;
extern Char_t S_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_S(Char_t*, Double_t, Double_t);
void Sort_S(Int_t l, Int_t r);
Int_t GetEnergyBins_S(Int_t* bins=NULL);
Int_t GetNPts_S();
Int_t ReadLine_S(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_S();
Double_t GetScale_S();

//-----------------------------------------------------------------------------

#endif
