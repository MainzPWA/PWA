#ifndef __Provide_Tz__
#define __Provide_Tz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Tz_val[EBINS][THBINS];
extern Double_t Tz_err[EBINS][THBINS];
extern Double_t Tz_unc[EBINS][THBINS];
extern Double_t Tz_th[EBINS][THBINS];
extern Double_t Tz_lo[EBINS];
extern Double_t Tz_en[EBINS];
extern Double_t Tz_hi[EBINS];
extern Double_t Tz_wt[EBINS];
extern Double_t Tz_sy[EBINS];
extern Double_t Tz_sc[EBINS];
extern Bool_t Tz_pre[EBINS];
extern Int_t Tz_pts[EBINS];
extern Int_t Tz_bin;
extern Char_t Tz_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_Tz(Char_t*, Double_t, Double_t);
void Sort_Tz(Int_t l, Int_t r);
Int_t GetEnergyBins_Tz(Int_t* bins=NULL);
Int_t GetNPts_Tz();
Int_t ReadLine_Tz(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Tz();
Double_t GetScale_Tz();

//-----------------------------------------------------------------------------

#endif
