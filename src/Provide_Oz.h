#ifndef __Provide_Oz__
#define __Provide_Oz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Oz_val[EBINS][THBINS];
extern Double_t Oz_err[EBINS][THBINS];
extern Double_t Oz_unc[EBINS][THBINS];
extern Double_t Oz_th[EBINS][THBINS];
extern Double_t Oz_lo[EBINS];
extern Double_t Oz_en[EBINS];
extern Double_t Oz_hi[EBINS];
extern Double_t Oz_wt[EBINS];
extern Double_t Oz_sy[EBINS];
extern Double_t Oz_sc[EBINS];
extern Bool_t Oz_pre[EBINS];
extern Int_t Oz_pts[EBINS];
extern Int_t Oz_bin;
extern Char_t Oz_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_Oz(Char_t*, Double_t, Double_t);
void Sort_Oz(Int_t l, Int_t r);
Int_t GetEnergyBins_Oz(Int_t* bins=NULL);
Int_t GetNPts_Oz();
Int_t ReadLine_Oz(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Oz();
Double_t GetScale_Oz();

//-----------------------------------------------------------------------------

#endif
