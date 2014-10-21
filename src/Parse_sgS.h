#ifndef __Parse_sgS__
#define __Parse_sgS__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgS_val[EBINS][THBINS];
extern Double_t sgS_err[EBINS][THBINS];
extern Double_t sgS_th[EBINS][THBINS];
extern Double_t sgS_en[EBINS];
extern Double_t sgS_wt[EBINS];
extern Double_t sgS_sy[EBINS];
extern Int_t sgS_pts[EBINS];
extern Int_t sgS_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgS();
Int_t GetEnergyBin_sgS();
Int_t ExistEnergyBin_sgS(Double_t);
Int_t ReadLine_sgS(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgS();
Double_t GetScale_sgS();

//-----------------------------------------------------------------------------

#endif

