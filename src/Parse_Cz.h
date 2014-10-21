#ifndef __Parse_Cz__
#define __Parse_Cz__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Cz_val[EBINS][THBINS];
extern Double_t Cz_err[EBINS][THBINS];
extern Double_t Cz_th[EBINS][THBINS];
extern Double_t Cz_en[EBINS];
extern Double_t Cz_wt[EBINS];
extern Double_t Cz_sy[EBINS];
extern Int_t Cz_pts[EBINS];
extern Int_t Cz_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_Cz();
Int_t GetEnergyBin_Cz();
Int_t ExistEnergyBin_Cz(Double_t);
Int_t ReadLine_Cz(FILE*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Cz();
Double_t GetScale_Cz();

//-----------------------------------------------------------------------------

#endif
