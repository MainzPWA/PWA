#ifndef __Parse_S__
#define __Parse_S__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t S_val[EBINS][THBINS];
extern Double_t S_err[EBINS][THBINS];
extern Double_t S_th[EBINS][THBINS];
extern Double_t S_en[EBINS];
extern Double_t S_wt[EBINS];
extern Double_t S_sy[EBINS];
extern Int_t S_pts[EBINS];
extern Int_t S_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_S();
Int_t GetEnergyBin_S();
Int_t ExistEnergyBin_S(Double_t);
Int_t ReadLine_S(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_S();
Double_t GetScale_S();

//-----------------------------------------------------------------------------

#endif

