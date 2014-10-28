#ifndef __Parse_sgH__
#define __Parse_sgH__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgH_val[EBINS][THBINS];
extern Double_t sgH_err[EBINS][THBINS];
extern Double_t sgH_th[EBINS][THBINS];
extern Double_t sgH_en[EBINS];
extern Double_t sgH_wt[EBINS];
extern Double_t sgH_sy[EBINS];
extern Int_t sgH_pts[EBINS];
extern Int_t sgH_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgH();
Int_t GetEnergyBin_sgH();
Int_t ExistEnergyBin_sgH(Double_t);
Int_t ReadLine_sgH(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgH();
Double_t GetScale_sgH();

//-----------------------------------------------------------------------------

#endif

