#ifndef __Provide_T__
#define __Provide_T__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t T_val[EBINS][THBINS];
extern Double_t T_err[EBINS][THBINS];
extern Double_t T_unc[EBINS][THBINS];
extern Double_t T_th[EBINS][THBINS];
extern Double_t T_lo[EBINS];
extern Double_t T_en[EBINS];
extern Double_t T_hi[EBINS];
extern Double_t T_wt[EBINS];
extern Double_t T_sy[EBINS];
extern Double_t T_sc[EBINS];
extern Bool_t T_pre[EBINS];
extern Int_t T_pts[EBINS];
extern Int_t T_bin;
extern Char_t T_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_T(Char_t*, Double_t, Double_t);
void Sort_T(Int_t l, Int_t r);
Int_t GetEnergyBins_T(Int_t* bins=NULL);
Int_t GetNPts_T();
Int_t ReadLine_T(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_T();
Double_t GetScale_T();

//-----------------------------------------------------------------------------

#endif
