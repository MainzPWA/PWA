#ifndef __Parse_E__
#define __Parse_E__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t E_val[EBINS][THBINS];
extern Double_t E_err[EBINS][THBINS];
extern Double_t E_th[EBINS][THBINS];
extern Double_t E_en[EBINS];
extern Double_t E_wt[EBINS];
extern Double_t E_sy[EBINS];
extern Int_t E_pts[EBINS];
extern Int_t E_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_E();
Int_t GetEnergyBin_E();
Int_t ExistEnergyBin_E(Double_t);
Int_t ReadLine_E(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_E();
Double_t GetScale_E();

//-----------------------------------------------------------------------------

#endif

