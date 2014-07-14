#ifndef __Parse_G__
#define __Parse_G__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t G_val[EBINS][THBINS];
extern Double_t G_err[EBINS][THBINS];
extern Double_t G_th[EBINS][THBINS];
extern Double_t G_en[EBINS];
extern Double_t G_wt[EBINS];
extern Double_t G_sy[EBINS];
extern Int_t G_pts[EBINS];
extern Int_t G_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_G();
Int_t GetEnergyBin_G();
Int_t ExistEnergyBin_G(Double_t);
Double_t GetChiSq_G();
Double_t GetScale_G();

//-----------------------------------------------------------------------------

#endif

