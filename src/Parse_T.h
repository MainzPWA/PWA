#ifndef __Parse_T__
#define __Parse_T__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t T_val[EBINS][THBINS];
extern Double_t T_err[EBINS][THBINS];
extern Double_t T_th[EBINS][THBINS];
extern Double_t T_en[EBINS];
extern Double_t T_wt[EBINS];
extern Double_t T_sy[EBINS];
extern Int_t T_pts[EBINS];
extern Int_t T_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_T();
Int_t GetEnergyBin_T();
Int_t ExistEnergyBin_T(Double_t);
Int_t ReadLine_T(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_T();
Double_t GetScale_T();

//-----------------------------------------------------------------------------

#endif

