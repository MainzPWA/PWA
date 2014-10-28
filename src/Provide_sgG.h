#ifndef __Provide_sgG__
#define __Provide_sgG__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgG_val[EBINS][THBINS];
extern Double_t sgG_err[EBINS][THBINS];
extern Double_t sgG_unc[EBINS][THBINS];
extern Double_t sgG_th[EBINS][THBINS];
extern Double_t sgG_lo[EBINS];
extern Double_t sgG_en[EBINS];
extern Double_t sgG_hi[EBINS];
extern Double_t sgG_wt[EBINS];
extern Double_t sgG_sy[EBINS];
extern Double_t sgG_sc[EBINS];
extern Bool_t sgG_pre[EBINS];
extern Int_t sgG_pts[EBINS];
extern Int_t sgG_bin;
extern Char_t sgG_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sgG(Char_t*, Double_t, Double_t);
void Sort_sgG(Int_t l, Int_t r);
Int_t GetEnergyBins_sgG(Int_t* bins=NULL);
Int_t GetNPts_sgG();
Int_t ReadLine_sgG(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgG();
Double_t GetScale_sgG();

//-----------------------------------------------------------------------------

#endif
