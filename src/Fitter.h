#ifndef __Fitter__
#define __Fitter__

#include "Constants.h"
#include "Globals.h"
#include "Provide_sg0.h"
#include "Provide_sgS.h"
#include "Provide_sgT.h"
#include "Provide_sgP.h"
#include "Provide_sgE.h"
#include "Provide_sgF.h"
#include "Provide_sgG.h"
#include "Provide_sgH.h"
#include "Provide_sgCx.h"
#include "Provide_sgCz.h"
#include "Provide_sgOx.h"
#include "Provide_sgOz.h"
#include "Provide_sgLx.h"
#include "Provide_sgLz.h"
#include "Provide_sgTx.h"
#include "Provide_sgTz.h"
#include "Provide_S.h"
#include "Provide_T.h"
#include "Provide_P.h"
#include "Provide_E.h"
#include "Provide_F.h"
#include "Provide_G.h"
#include "Provide_H.h"
#include "Provide_Cx.h"
#include "Provide_Cz.h"
#include "Provide_Ox.h"
#include "Provide_Oz.h"
#include "Provide_Lx.h"
#include "Provide_Lz.h"
#include "Provide_Tx.h"
#include "Provide_Tz.h"
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
Double_t PenaltyHELI();
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
    //Sum up squared model amplitude magnitudes (will be used as normalisation factor)
    MagSq+=F_maid.Rho2();
  }
  return SumSq/MagSq;
}

//-----------------------------------------------------------------------------

inline Double_t EvaluateHELI(Double_t CosTheta, Int_t eM)
{
  Double_t SumSq = 0.0;
  Double_t MagSq = 0.0;
  TComplex H_data;
  TComplex H_maid;

  //Sum up H1...H4 amplitudes
  for(Int_t i=1; i<5; i++)
  {
    //Get amplitudes from current fit parameters and model values
    H_data = H(i, CosTheta);
    H_maid = maid_H(i, CosTheta, eM);
    //Sum up all squared amplitude deviations with adjustable weight factor
    SumSq+=(H_data - H_maid).Rho2()/(WEIGHT[i]*WEIGHT[i]);
    //Sum up squared model amplitude magnitudes (will be used as normalisation factor)
    MagSq+=H_maid.Rho2();
  }
  return SumSq/MagSq;
}

//-----------------------------------------------------------------------------

#endif

