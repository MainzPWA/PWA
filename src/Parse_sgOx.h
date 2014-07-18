#ifndef __Parse_sgOx__
#define __Parse_sgOx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgOx_val[EBINS][THBINS];
extern Double_t sgOx_err[EBINS][THBINS];
extern Double_t sgOx_th[EBINS][THBINS];
extern Double_t sgOx_en[EBINS];
extern Double_t sgOx_wt[EBINS];
extern Double_t sgOx_sy[EBINS];
extern Int_t sgOx_pts[EBINS];
extern Int_t sgOx_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgOx();
Int_t GetEnergyBin_sgOx();
Int_t ExistEnergyBin_sgOx(Double_t);
Double_t GetChiSq_sgOx();
Double_t GetScale_sgOx();

//-----------------------------------------------------------------------------

#endif
