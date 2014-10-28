#ifndef __Provide_sgH__
#define __Provide_sgH__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgH_val[EBINS][THBINS];
extern Double_t sgH_err[EBINS][THBINS];
extern Double_t sgH_unc[EBINS][THBINS];
extern Double_t sgH_th[EBINS][THBINS];
extern Double_t sgH_lo[EBINS];
extern Double_t sgH_en[EBINS];
extern Double_t sgH_hi[EBINS];
extern Double_t sgH_wt[EBINS];
extern Double_t sgH_sy[EBINS];
extern Double_t sgH_sc[EBINS];
extern Bool_t sgH_pre[EBINS];
extern Int_t sgH_pts[EBINS];
extern Int_t sgH_bin;
extern Char_t sgH_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgH(Char_t*, Double_t, Double_t);
void Sort_sgH(Int_t l, Int_t r);
Int_t GetEnergyBins_sgH(Int_t* bins=NULL);
Int_t GetNPts_sgH();
Int_t ReadLine_sgH(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgH();
Double_t GetScale_sgH();

//-----------------------------------------------------------------------------

#endif
