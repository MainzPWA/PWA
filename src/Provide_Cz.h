#ifndef __Provide_Cz__
#define __Provide_Cz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Cz_val[EBINS][THBINS];
extern Double_t Cz_err[EBINS][THBINS];
extern Double_t Cz_unc[EBINS][THBINS];
extern Double_t Cz_th[EBINS][THBINS];
extern Double_t Cz_lo[EBINS];
extern Double_t Cz_en[EBINS];
extern Double_t Cz_hi[EBINS];
extern Double_t Cz_wt[EBINS];
extern Double_t Cz_sy[EBINS];
extern Double_t Cz_sc[EBINS];
extern Bool_t Cz_pre[EBINS];
extern Int_t Cz_pts[EBINS];
extern Int_t Cz_bin;
extern Char_t Cz_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_Cz(Char_t*, Double_t, Double_t);
void Sort_Cz(Int_t l, Int_t r);
Int_t GetEnergyBins_Cz(Int_t* bins=NULL);
Int_t GetNPts_Cz();
Int_t ReadLine_Cz(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Cz();
Double_t GetScale_Cz();

//-----------------------------------------------------------------------------

#endif
