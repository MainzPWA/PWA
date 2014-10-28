#ifndef __Provide_G__
#define __Provide_G__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t G_val[EBINS][THBINS];
extern Double_t G_err[EBINS][THBINS];
extern Double_t G_unc[EBINS][THBINS];
extern Double_t G_th[EBINS][THBINS];
extern Double_t G_lo[EBINS];
extern Double_t G_en[EBINS];
extern Double_t G_hi[EBINS];
extern Double_t G_wt[EBINS];
extern Double_t G_sy[EBINS];
extern Double_t G_sc[EBINS];
extern Bool_t G_pre[EBINS];
extern Int_t G_pts[EBINS];
extern Int_t G_bin;
extern Char_t G_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_G(Char_t*, Double_t, Double_t);
void Sort_G(Int_t l, Int_t r);
Int_t GetEnergyBins_G(Int_t* bins=NULL);
Int_t GetNPts_G();
Int_t ReadLine_G(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_G();
Double_t GetScale_G();

//-----------------------------------------------------------------------------

#endif
