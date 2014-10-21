#ifndef __Parse_sgG__
#define __Parse_sgG__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sgG_val[EBINS][THBINS];
extern Double_t sgG_err[EBINS][THBINS];
extern Double_t sgG_th[EBINS][THBINS];
extern Double_t sgG_en[EBINS];
extern Double_t sgG_wt[EBINS];
extern Double_t sgG_sy[EBINS];
extern Int_t sgG_pts[EBINS];
extern Int_t sgG_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sgG();
Int_t GetEnergyBin_sgG();
Int_t ExistEnergyBin_sgG(Double_t);
Int_t ReadLine_sgG(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sgG();
Double_t GetScale_sgG();

//-----------------------------------------------------------------------------

#endif

