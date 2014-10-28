#ifndef __Provide_Ox__
#define __Provide_Ox__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Ox_val[EBINS][THBINS];
extern Double_t Ox_err[EBINS][THBINS];
extern Double_t Ox_unc[EBINS][THBINS];
extern Double_t Ox_th[EBINS][THBINS];
extern Double_t Ox_lo[EBINS];
extern Double_t Ox_en[EBINS];
extern Double_t Ox_hi[EBINS];
extern Double_t Ox_wt[EBINS];
extern Double_t Ox_sy[EBINS];
extern Double_t Ox_sc[EBINS];
extern Bool_t Ox_pre[EBINS];
extern Int_t Ox_pts[EBINS];
extern Int_t Ox_bin;
extern Char_t Ox_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_Ox(Char_t*, Double_t, Double_t);
void Sort_Ox(Int_t l, Int_t r);
Int_t GetEnergyBins_Ox(Int_t* bins=NULL);
Int_t GetNPts_Ox();
Int_t ReadLine_Ox(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Ox();
Double_t GetScale_Ox();

//-----------------------------------------------------------------------------

#endif
