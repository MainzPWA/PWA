#ifndef __Parse_F__
#define __Parse_F__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t F_val[EBINS][THBINS];
extern Double_t F_err[EBINS][THBINS];
extern Double_t F_th[EBINS][THBINS];
extern Double_t F_en[EBINS];
extern Double_t F_wt[EBINS];
extern Double_t F_sy[EBINS];
extern Int_t F_pts[EBINS];
extern Int_t F_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_F();
Int_t GetEnergyBin_F();
Int_t ExistEnergyBin_F(Double_t);
Double_t GetChiSq_F();
Double_t GetScale_F();

//-----------------------------------------------------------------------------

#endif

