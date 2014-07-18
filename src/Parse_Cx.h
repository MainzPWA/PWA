#ifndef __Parse_Cx__
#define __Parse_Cx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Cx_val[EBINS][THBINS];
extern Double_t Cx_err[EBINS][THBINS];
extern Double_t Cx_th[EBINS][THBINS];
extern Double_t Cx_en[EBINS];
extern Double_t Cx_wt[EBINS];
extern Double_t Cx_sy[EBINS];
extern Int_t Cx_pts[EBINS];
extern Int_t Cx_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_Cx();
Int_t GetEnergyBin_Cx();
Int_t ExistEnergyBin_Cx(Double_t);
Double_t GetChiSq_Cx();
Double_t GetScale_Cx();

//-----------------------------------------------------------------------------

#endif
