#ifndef __Provide_F__
#define __Provide_F__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t F_val[EBINS][THBINS];
extern Double_t F_err[EBINS][THBINS];
extern Double_t F_unc[EBINS][THBINS];
extern Double_t F_th[EBINS][THBINS];
extern Double_t F_lo[EBINS];
extern Double_t F_en[EBINS];
extern Double_t F_hi[EBINS];
extern Double_t F_wt[EBINS];
extern Double_t F_sy[EBINS];
extern Double_t F_sc[EBINS];
extern Bool_t F_pre[EBINS];
extern Int_t F_pts[EBINS];
extern Int_t F_bin;
extern Char_t F_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_F(Char_t*, Double_t, Double_t);
void Sort_F(Int_t l, Int_t r);
Int_t GetEnergyBins_F(Int_t* bins=NULL);
Int_t GetNPts_F();
Int_t ReadLine_F(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_F();
Double_t GetScale_F();

//-----------------------------------------------------------------------------

#endif
