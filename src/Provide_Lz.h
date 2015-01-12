#ifndef __Provide_Lz__
#define __Provide_Lz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Lz_val[EBINS][THBINS];
extern Double_t Lz_err[EBINS][THBINS];
extern Double_t Lz_unc[EBINS][THBINS];
extern Double_t Lz_th[EBINS][THBINS];
extern Double_t Lz_lo[EBINS];
extern Double_t Lz_en[EBINS];
extern Double_t Lz_hi[EBINS];
extern Double_t Lz_wt[EBINS];
extern Double_t Lz_sy[EBINS];
extern Double_t Lz_sc[EBINS];
extern Bool_t Lz_pre[EBINS];
extern Int_t Lz_pts[EBINS];
extern Int_t Lz_bin;
extern Char_t Lz_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_Lz(Char_t*, Double_t, Double_t);
void Sort_Lz(Int_t l, Int_t r);
Int_t GetEnergyBins_Lz(Int_t* bins=NULL);
Int_t GetNPts_Lz();
Int_t ReadLine_Lz(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Lz();
Double_t GetScale_Lz();

//-----------------------------------------------------------------------------

#endif
