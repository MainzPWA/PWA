#ifndef __Parse_sgCx__
#define __Parse_sgCx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgCx_val[EBINS][THBINS];
extern Double_t sgCx_err[EBINS][THBINS];
extern Double_t sgCx_th[EBINS][THBINS];
extern Double_t sgCx_en[EBINS];
extern Double_t sgCx_wt[EBINS];
extern Double_t sgCx_sy[EBINS];
extern Int_t sgCx_pts[EBINS];
extern Int_t sgCx_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgCx();
Int_t GetEnergyBin_sgCx();
Int_t ExistEnergyBin_sgCx(Double_t);
Int_t ReadLine_sgCx(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgCx();
Double_t GetScale_sgCx();

//-----------------------------------------------------------------------------

#endif
