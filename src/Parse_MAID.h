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

//MAID CGLN amplitude F1 in expansion up to l_max
inline TComplex maid_F1(Double_t CosTheta, Int_t e)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=0; l<L_MAX+1; l++)
    Complex+=((maid_Mp[l][e]*(1.0*l) + maid_Ep[l][e]) * DL(l+1, CosTheta) + (maid_Mm[l][e]*(1.0*l + 1.0) + maid_Em[l][e]) * DL(l-1, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//MAID CGLN amplitude F2 in expansion up to l_max
inline TComplex maid_F2(Double_t CosTheta, Int_t e)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=1; l<L_MAX+1; l++)
    Complex+=((maid_Mp[l][e]*(1.0*l + 1.0) + maid_Mm[l][e]*(1.0*l)) * DL(l, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//MAID CGLN amplitude F3 in expansion up to l_max
inline TComplex maid_F3(Double_t CosTheta, Int_t e)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=1; l<L_MAX+1; l++)
    Complex+=((maid_Ep[l][e] - maid_Mp[l][e]) * D2L(l+1, CosTheta) + (maid_Em[l][e] + maid_Mm[l][e]) * D2L(l-1, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//MAID CGLN amplitude F4 in expansion up to l_max
inline TComplex maid_F4(Double_t CosTheta, Int_t e)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=2; l<L_MAX+1; l++)
    Complex+=((maid_Mp[l][e] - maid_Ep[l][e] - maid_Mm[l][e] - maid_Em[l][e]) * D2L(l, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//Wrapper for explicit MAID CGLN amplitudes F_i
inline TComplex maid_F(Int_t i, Double_t CosTheta, Int_t e)
{
  switch(i)
  {
    case 1: return maid_F1(CosTheta, e);
    case 2: return maid_F2(CosTheta, e);
    case 3: return maid_F3(CosTheta, e);
    case 4: return maid_F4(CosTheta, e);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

inline TComplex maid_F1cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_F1(CosTheta, e)); }
inline TComplex maid_F2cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_F2(CosTheta, e)); }
inline TComplex maid_F3cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_F3(CosTheta, e)); }
inline TComplex maid_F4cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(maid_F4(CosTheta, e)); }

//-----------------------------------------------------------------------------

//Wrapper for explicit MAID CGLN amplitudes F_i*
inline TComplex maid_Fcc(Int_t i, Double_t CosTheta, Int_t e)
{
  switch(i)
  {
    case 1: return maid_F1cc(CosTheta, e);
    case 2: return maid_F2cc(CosTheta, e);
    case 3: return maid_F3cc(CosTheta, e);
    case 4: return maid_F4cc(CosTheta, e);
    default: return TComplex(0.0, 0.0);
  }
}

//-----------------------------------------------------------------------------

#endif

