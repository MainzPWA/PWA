#ifndef __Fitter__
#define __Fitter__

#include "Constants.h"
#include "Globals.h"
#include "Parse_sg0.h"
#include "Parse_sgS.h"
#include "Parse_sgT.h"
#include "Parse_sgP.h"
#include "Parse_sgE.h"
#include "Parse_sgF.h"
#include "Parse_sgG.h"
#include "Parse_sgH.h"
#include "Parse_sgCx.h"
#include "Parse_sgCz.h"
#include "Parse_sgOx.h"
#include "Parse_sgOz.h"
#include "Parse_S.h"
#include "Parse_T.h"
#include "Parse_P.h"
#include "Parse_E.h"
#include "Parse_F.h"
#include "Parse_G.h"
#include "Parse_H.h"
#include "Parse_Cx.h"
#include "Parse_Cz.h"
#include "Parse_Ox.h"
#include "Parse_Oz.h"
#include "Parse_MAID.h"
#include "Multipoles.h"

//-----------------------------------------------------------------------------

extern TMinuit* gMinuit;

//Errors on multipoles
extern TComplex DEp[LBINS];
extern TComplex DEm[LBINS];
extern TComplex DMp[LBINS];
extern TComplex DMm[LBINS];

extern Double_t f_obs[OBS];
extern Double_t Df_obs[OBS];

//-----------------------------------------------------------------------------

void fcn_main(Double_t *par);
void fcn_chi_pen(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcn_chi(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void InitMinuit();
void Fit();
void SortFits(Double_t* Quality, Double_t* MeanErr, Int_t* Attempt, Int_t l, Int_t r);
void SetSPWaves();
void SetWaves(Int_t l);
void SetRealPWaves();
void SetUnitarity();
void FixSPWaves();
void FixSPPhases();
void FixPhases(Int_t l);
void FixWaves(Int_t l);
void FixRealPWaves();
void FixThreshold();
void FixUnitarity();
void FixScalings();
void SetParameters();
void GetParameters(Double_t* Par, Double_t* Err);
void UseParameters(Double_t* Par, Double_t* Err);
void StoreFit(Int_t);
void StoreModel();
Bool_t SingleFit();
Double_t ChiSq();
Double_t Scale();
Double_t Penalty();
Double_t PenaltyMLP1();
Double_t PenaltyMLP2();
Double_t PenaltyCGLN();
Double_t VariateRel();
Double_t VariateAbs();
Double_t GetErrors(Double_t* Par, Double_t* Err);
Double_t ErrorBase(Double_t* Par, Double_t* Err);
Int_t NPts();
Int_t NPar();
Int_t NMlp();
Int_t NSca();
Int_t NDF();

//-----------------------------------------------------------------------------

inline Double_t ErrorChi2(Double_t* Par, Double_t* Err)
{
  //Initialise Minuit fitter
  InitMinuit();
  //Use fcn without penalty
  gMinuit->SetFCN(fcn_chi);
  return ErrorBase(Par, Err);
}

//-----------------------------------------------------------------------------

inline Double_t ErrorChi2Penalty(Double_t* Par, Double_t* Err)
{
  //Initialise Minuit fitter
  InitMinuit();
  //Use fcn with penalty
  gMinuit->SetFCN(fcn_chi_pen);
  return ErrorBase(Par, Err);
}

//-----------------------------------------------------------------------------

inline Double_t EvaluateCGLN(Double_t CosTheta, Int_t eM)
{
  Double_t SumSq = 0.0;
  Double_t MagSq = 0.0;
  TComplex F_data;
  TComplex F_maid;

  //Sum up F1...F4 amplitudes
  for(Int_t i=1; i<5; i++)
  {
    //Get amplitudes from current fit parameters and model values
    F_data = F(i, CosTheta);
    F_maid = maid_F(i, CosTheta, eM);
    //Sum up all squared amplitude deviations with adjustable weight factor
    SumSq+=(F_data - F_maid).Rho2()/(WEIGHT[i]*WEIGHT[i]);
    //Sum up model amplitude magnitudes (will be used as normalisation factor)
    MagSq+=F_maid.Rho2();
  }
  return SumSq/MagSq;
}

//-----------------------------------------------------------------------------

#endif

