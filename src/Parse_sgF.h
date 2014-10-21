#ifndef __Parse_sgF__
#define __Parse_sgF__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgF_val[EBINS][THBINS];
extern Double_t sgF_err[EBINS][THBINS];
extern Double_t sgF_th[EBINS][THBINS];
extern Double_t sgF_en[EBINS];
extern Double_t sgF_wt[EBINS];
extern Double_t sgF_sy[EBINS];
extern Int_t sgF_pts[EBINS];
extern Int_t sgF_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgF();
Int_t GetEnergyBin_sgF();
Int_t ExistEnergyBin_sgF(Double_t);
Int_t ReadLine_sgF(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgF();
Double_t GetScale_sgF();

//-----------------------------------------------------------------------------

#endif
