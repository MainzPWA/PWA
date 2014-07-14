#ifndef __Parse_sgP__
#define __Parse_sgP__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgP_val[EBINS][THBINS];
extern Double_t sgP_err[EBINS][THBINS];
extern Double_t sgP_th[EBINS][THBINS];
extern Double_t sgP_en[EBINS];
extern Double_t sgP_wt[EBINS];
extern Double_t sgP_sy[EBINS];
extern Int_t sgP_pts[EBINS];
extern Int_t sgP_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgP();
Int_t GetEnergyBin_sgP();
Int_t ExistEnergyBin_sgP(Double_t);
Double_t GetChiSq_sgP();
Double_t GetScale_sgP();

//-----------------------------------------------------------------------------

#endif

