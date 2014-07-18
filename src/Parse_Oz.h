#ifndef __Parse_Oz__
#define __Parse_Oz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Oz_val[EBINS][THBINS];
extern Double_t Oz_err[EBINS][THBINS];
extern Double_t Oz_th[EBINS][THBINS];
extern Double_t Oz_en[EBINS];
extern Double_t Oz_wt[EBINS];
extern Double_t Oz_sy[EBINS];
extern Int_t Oz_pts[EBINS];
extern Int_t Oz_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_Oz();
Int_t GetEnergyBin_Oz();
Int_t ExistEnergyBin_Oz(Double_t);
Double_t GetChiSq_Oz();
Double_t GetScale_Oz();

//-----------------------------------------------------------------------------

#endif
