#ifndef __Provide_sg0__
#define __Provide_sg0__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t sg0_val[EBINS][THBINS];
extern Double_t sg0_err[EBINS][THBINS];
extern Double_t sg0_unc[EBINS][THBINS];
extern Double_t sg0_th[EBINS][THBINS];
extern Double_t sg0_lo[EBINS];
extern Double_t sg0_en[EBINS];
extern Double_t sg0_hi[EBINS];
extern Double_t sg0_wt[EBINS];
extern Double_t sg0_sy[EBINS];
extern Double_t sg0_sc[EBINS];
extern Bool_t sg0_pre[EBINS];
extern Int_t sg0_pts[EBINS];
extern Int_t sg0_bin;
extern Char_t sg0_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_sg0(Char_t*, Double_t, Double_t);
void Sort_sg0(Int_t l, Int_t r);
Int_t GetEnergyBins_sg0(Int_t* bins=NULL);
Int_t GetNPts_sg0();
Int_t ReadLine_sg0(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_sg0();
Double_t GetScale_sg0();

//-----------------------------------------------------------------------------

#endif
