#ifndef __Provide_Cx__
#define __Provide_Cx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Cx_val[EBINS][THBINS];
extern Double_t Cx_err[EBINS][THBINS];
extern Double_t Cx_unc[EBINS][THBINS];
extern Double_t Cx_th[EBINS][THBINS];
extern Double_t Cx_lo[EBINS];
extern Double_t Cx_en[EBINS];
extern Double_t Cx_hi[EBINS];
extern Double_t Cx_wt[EBINS];
extern Double_t Cx_sy[EBINS];
extern Double_t Cx_sc[EBINS];
extern Bool_t Cx_pre[EBINS];
extern Int_t Cx_pts[EBINS];
extern Int_t Cx_bin;
extern Char_t Cx_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_Cx(Char_t*, Double_t, Double_t);
void Sort_Cx(Int_t l, Int_t r);
Int_t GetEnergyBins_Cx(Int_t* bins=NULL);
Int_t GetNPts_Cx();
Int_t ReadLine_Cx(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Cx();
Double_t GetScale_Cx();

//-----------------------------------------------------------------------------

#endif
