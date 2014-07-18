#ifndef __Parse_Ox__
#define __Parse_Ox__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Ox_val[EBINS][THBINS];
extern Double_t Ox_err[EBINS][THBINS];
extern Double_t Ox_th[EBINS][THBINS];
extern Double_t Ox_en[EBINS];
extern Double_t Ox_wt[EBINS];
extern Double_t Ox_sy[EBINS];
extern Int_t Ox_pts[EBINS];
extern Int_t Ox_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_Ox();
Int_t GetEnergyBin_Ox();
Int_t ExistEnergyBin_Ox(Double_t);
Double_t GetChiSq_Ox();
Double_t GetScale_Ox();

//-----------------------------------------------------------------------------

#endif
