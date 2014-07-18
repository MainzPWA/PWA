#ifndef __Parse_sgCz__
#define __Parse_sgCz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgCz_val[EBINS][THBINS];
extern Double_t sgCz_err[EBINS][THBINS];
extern Double_t sgCz_th[EBINS][THBINS];
extern Double_t sgCz_en[EBINS];
extern Double_t sgCz_wt[EBINS];
extern Double_t sgCz_sy[EBINS];
extern Int_t sgCz_pts[EBINS];
extern Int_t sgCz_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgCz();
Int_t GetEnergyBin_sgCz();
Int_t ExistEnergyBin_sgCz(Double_t);
Double_t GetChiSq_sgCz();
Double_t GetScale_sgCz();

//-----------------------------------------------------------------------------

#endif
