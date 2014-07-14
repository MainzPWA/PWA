#ifndef __Parse_MAID__
#define __Parse_MAID__

#include "Constants.h"
#include "Globals.h"
#include "Legendre.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern TComplex maid_Ep[LBINS][EBINS];
extern TComplex maid_Em[LBINS][EBINS];
extern TComplex maid_Mp[LBINS][EBINS];
extern TComplex maid_Mm[LBINS][EBINS];
extern Double_t maid_en[EBINS];
extern Int_t maid_bin;

//-----------------------------------------------------------------------------

void Parse_MAID();
Int_t GetEnergyBin_maid();

//-----------------------------------------------------------------------------

#endif

