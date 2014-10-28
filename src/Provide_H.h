#ifndef __Provide_H__
#define __Provide_H__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t H_val[EBINS][THBINS];
extern Double_t H_err[EBINS][THBINS];
extern Double_t H_unc[EBINS][THBINS];
extern Double_t H_th[EBINS][THBINS];
extern Double_t H_lo[EBINS];
extern Double_t H_en[EBINS];
extern Double_t H_hi[EBINS];
extern Double_t H_wt[EBINS];
extern Double_t H_sy[EBINS];
extern Double_t H_sc[EBINS];
extern Bool_t H_pre[EBINS];
extern Int_t H_pts[EBINS];
extern Int_t H_bin;
extern Char_t H_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_H(Char_t*, Double_t, Double_t);
void Sort_H(Int_t l, Int_t r);
Int_t GetEnergyBins_H(Int_t* bins=NULL);
Int_t GetNPts_H();
Int_t ReadLine_H(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_H();
Double_t GetScale_H();

//-----------------------------------------------------------------------------

#endif
