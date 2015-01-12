#ifndef __Provide_Lx__
#define __Provide_Lx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Lx_val[EBINS][THBINS];
extern Double_t Lx_err[EBINS][THBINS];
extern Double_t Lx_unc[EBINS][THBINS];
extern Double_t Lx_th[EBINS][THBINS];
extern Double_t Lx_lo[EBINS];
extern Double_t Lx_en[EBINS];
extern Double_t Lx_hi[EBINS];
extern Double_t Lx_wt[EBINS];
extern Double_t Lx_sy[EBINS];
extern Double_t Lx_sc[EBINS];
extern Bool_t Lx_pre[EBINS];
extern Int_t Lx_pts[EBINS];
extern Int_t Lx_bin;
extern Char_t Lx_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_Lx(Char_t*, Double_t, Double_t);
void Sort_Lx(Int_t l, Int_t r);
Int_t GetEnergyBins_Lx(Int_t* bins=NULL);
Int_t GetNPts_Lx();
Int_t ReadLine_Lx(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Lx();
Double_t GetScale_Lx();

//-----------------------------------------------------------------------------

#endif
