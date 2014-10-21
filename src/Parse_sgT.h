#ifndef __Parse_sgT__
#define __Parse_sgT__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgT_val[EBINS][THBINS];
extern Double_t sgT_err[EBINS][THBINS];
extern Double_t sgT_th[EBINS][THBINS];
extern Double_t sgT_en[EBINS];
extern Double_t sgT_wt[EBINS];
extern Double_t sgT_sy[EBINS];
extern Int_t sgT_pts[EBINS];
extern Int_t sgT_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgT();
void Sort_sgT(Int_t l, Int_t r);
Int_t GetEnergyBin_sgT();
Int_t ExistEnergyBin_sgT(Double_t);
Int_t ReadLine_sgT(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgT();
Double_t GetScale_sgT();

//-----------------------------------------------------------------------------

#endif
