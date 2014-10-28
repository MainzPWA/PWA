#ifndef __Provide_E__
#define __Provide_E__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t E_val[EBINS][THBINS];
extern Double_t E_err[EBINS][THBINS];
extern Double_t E_unc[EBINS][THBINS];
extern Double_t E_th[EBINS][THBINS];
extern Double_t E_lo[EBINS];
extern Double_t E_en[EBINS];
extern Double_t E_hi[EBINS];
extern Double_t E_wt[EBINS];
extern Double_t E_sy[EBINS];
extern Double_t E_sc[EBINS];
extern Bool_t E_pre[EBINS];
extern Int_t E_pts[EBINS];
extern Int_t E_bin;
extern Char_t E_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_E(Char_t*, Double_t, Double_t);
void Sort_E(Int_t l, Int_t r);
Int_t GetEnergyBins_E(Int_t* bins=NULL);
Int_t GetNPts_E();
Int_t ReadLine_E(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_E();
Double_t GetScale_E();

//-----------------------------------------------------------------------------

#endif
