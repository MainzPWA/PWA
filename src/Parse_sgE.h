#ifndef __Parse_sgE__
#define __Parse_sgE__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgE_val[EBINS][THBINS];
extern Double_t sgE_err[EBINS][THBINS];
extern Double_t sgE_th[EBINS][THBINS];
extern Double_t sgE_en[EBINS];
extern Double_t sgE_wt[EBINS];
extern Double_t sgE_sy[EBINS];
extern Int_t sgE_pts[EBINS];
extern Int_t sgE_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgE();
Int_t GetEnergyBin_sgE();
Int_t ExistEnergyBin_sgE(Double_t);
Int_t ReadLine_sgE(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgE();
Double_t GetScale_sgE();

//-----------------------------------------------------------------------------

#endif

