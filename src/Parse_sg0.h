#ifndef __Parse_sg0__
#define __Parse_sg0__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sg0_val[EBINS][THBINS];
extern Double_t sg0_err[EBINS][THBINS];
extern Double_t sg0_th[EBINS][THBINS];
extern Double_t sg0_en[EBINS];
extern Double_t sg0_wt[EBINS];
extern Double_t sg0_sy[EBINS];
extern Int_t sg0_pts[EBINS];
extern Int_t sg0_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_sg0();
void Sort_sg0(Int_t l, Int_t r);
Int_t GetEnergyBin_sg0();
Int_t ExistEnergyBin_sg0(Double_t);
Int_t ReadLine_sg0(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sg0();
Double_t GetScale_sg0();

//-----------------------------------------------------------------------------

#endif
