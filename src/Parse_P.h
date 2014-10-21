#ifndef __Parse_P__
#define __Parse_P__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t P_val[EBINS][THBINS];
extern Double_t P_err[EBINS][THBINS];
extern Double_t P_th[EBINS][THBINS];
extern Double_t P_en[EBINS];
extern Double_t P_wt[EBINS];
extern Double_t P_sy[EBINS];
extern Int_t P_pts[EBINS];
extern Int_t P_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_P();
Int_t GetEnergyBin_P();
Int_t ExistEnergyBin_P(Double_t);
Int_t ReadLine_P(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_P();
Double_t GetScale_P();

//-----------------------------------------------------------------------------

#endif

