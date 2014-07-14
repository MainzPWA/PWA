#ifndef __Parse_H__
#define __Parse_H__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t H_val[EBINS][THBINS];
extern Double_t H_err[EBINS][THBINS];
extern Double_t H_th[EBINS][THBINS];
extern Double_t H_en[EBINS];
extern Double_t H_wt[EBINS];
extern Double_t H_sy[EBINS];
extern Int_t H_pts[EBINS];
extern Int_t H_bin;

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Parse_H();
Int_t GetEnergyBin_H();
Int_t ExistEnergyBin_H(Double_t);
Double_t GetChiSq_H();
Double_t GetScale_H();

//-----------------------------------------------------------------------------

#endif

