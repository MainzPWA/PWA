#ifndef __Provide_P__
#define __Provide_P__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t P_val[EBINS][THBINS];
extern Double_t P_err[EBINS][THBINS];
extern Double_t P_unc[EBINS][THBINS];
extern Double_t P_th[EBINS][THBINS];
extern Double_t P_lo[EBINS];
extern Double_t P_en[EBINS];
extern Double_t P_hi[EBINS];
extern Double_t P_wt[EBINS];
extern Double_t P_sy[EBINS];
extern Double_t P_sc[EBINS];
extern Bool_t P_pre[EBINS];
extern Int_t P_pts[EBINS];
extern Int_t P_bin;
extern Char_t P_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_P(Char_t*, Double_t, Double_t);
void Sort_P(Int_t l, Int_t r);
Int_t GetEnergyBins_P(Int_t* bins=NULL);
Int_t GetNPts_P();
Int_t ReadLine_P(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_P();
Double_t GetScale_P();

//-----------------------------------------------------------------------------

#endif
