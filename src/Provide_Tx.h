#ifndef __Provide_Tx__
#define __Provide_Tx__

#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern Double_t Tx_val[EBINS][THBINS];
extern Double_t Tx_err[EBINS][THBINS];
extern Double_t Tx_unc[EBINS][THBINS];
extern Double_t Tx_th[EBINS][THBINS];
extern Double_t Tx_lo[EBINS];
extern Double_t Tx_en[EBINS];
extern Double_t Tx_hi[EBINS];
extern Double_t Tx_wt[EBINS];
extern Double_t Tx_sy[EBINS];
extern Double_t Tx_sc[EBINS];
extern Bool_t Tx_pre[EBINS];
extern Int_t Tx_pts[EBINS];
extern Int_t Tx_bin;
extern Char_t Tx_id[EBINS][256];

extern Double_t f_obs[OBS];

//-----------------------------------------------------------------------------

void Load_Tx(Char_t*, Double_t, Double_t);
void Sort_Tx(Int_t l, Int_t r);
Int_t GetEnergyBins_Tx(Int_t* bins=NULL);
Int_t GetNPts_Tx();
Int_t ReadLine_Tx(FILE*, Double_t*, Double_t*, Double_t*, Double_t*);
Double_t GetChiSq_Tx();
Double_t GetScale_Tx();

//-----------------------------------------------------------------------------

#endif
