#ifndef __Parse_sgOz__
#define __Parse_sgOz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgOz_val[EBINS][THBINS];
extern Double_t sgOz_err[EBINS][THBINS];
extern Double_t sgOz_th[EBINS][THBINS];
extern Double_t sgOz_en[EBINS];
extern Double_t sgOz_wt[EBINS];
extern Double_t sgOz_sy[EBINS];
extern Int_t sgOz_pts[EBINS];
extern Int_t sgOz_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgOz();
Int_t GetEnergyBin_sgOz();
Int_t ExistEnergyBin_sgOz(Double_t);
Int_t ReadLine_sgOz(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgOz();
Double_t GetScale_sgOz();

//-----------------------------------------------------------------------------

#endif
