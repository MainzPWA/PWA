#ifndef __Constants__
#define __Constants__

//-----------------------------------------------------------------------------

#include "Root.h"

//-----------------------------------------------------------------------------

using namespace TMath;

//-----------------------------------------------------------------------------

static const Double_t MASS_PROTON      =  938.2720; //MeV
static const Double_t MASS_NEUTRON     =  939.5654; //MeV
static const Double_t MASS_LAMBDA      = 1115.6830; //MeV
static const Double_t MASS_SIGMAPLUS   = 1189.3700; //MeV
static const Double_t MASS_SIGMAZERO   = 1192.6420; //MeV
static const Double_t MASS_SIGMAMINUS  = 1197.4490; //MeV
static const Double_t MASS_PIZERO      =  134.9766; //MeV
static const Double_t MASS_PIPLUS      =  139.5702; //MeV
static const Double_t MASS_PIMINUS     =  139.5702; //MeV
static const Double_t MASS_KPLUS       =  493.6770; //MeV
static const Double_t MASS_KMINUS      =  493.6770; //MeV
static const Double_t MASS_KZERO       =  497.6140; //MeV
static const Double_t MASS_ETA         =  547.8530; //MeV
static const Double_t MASS_ETAPRIME    =  957.7800; //MeV
static const Double_t MASS2_PROTON     = MASS_PROTON*MASS_PROTON;
static const Double_t MASS2_NEUTRON    = MASS_NEUTRON*MASS_NEUTRON;
static const Double_t MASS2_LAMBDA     = MASS_LAMBDA*MASS_LAMBDA;
static const Double_t MASS2_SIGMAPLUS  = MASS_SIGMAPLUS*MASS_SIGMAPLUS;
static const Double_t MASS2_SIGMAZERO  = MASS_SIGMAZERO*MASS_SIGMAZERO;
static const Double_t MASS2_SIGMAMINUS = MASS_SIGMAMINUS*MASS_SIGMAMINUS;
static const Double_t MASS2_PIZERO     = MASS_PIZERO*MASS_PIZERO;
static const Double_t MASS2_PIMINUS    = MASS_PIMINUS*MASS_PIMINUS;
static const Double_t MASS2_PIPLUS     = MASS_PIPLUS*MASS_PIPLUS;
static const Double_t MASS2_KPLUS      = MASS_KPLUS*MASS_KPLUS;
static const Double_t MASS2_KMINUS     = MASS_KMINUS*MASS_KMINUS;
static const Double_t MASS2_KZERO      = MASS_KZERO*MASS_KZERO;
static const Double_t MASS2_ETA        = MASS_ETA*MASS_ETA;
static const Double_t MASS2_ETAPRIME   = MASS_ETAPRIME*MASS_ETAPRIME;
static const Double_t HBARC = 197.3269; //MeV fm (hbar*c)
static const Double_t UNIT  = (HBARC/MASS_PIPLUS)*(HBARC/MASS_PIPLUS)*1e-2;
//'Unit' is the conversion from the multipole's units (10^-3/m_pi+) to cross section unit (ub).
//Cross sections are proportional to the square of multipoles.

enum{ REL = 0, ABS = 1 };
enum{ MODEL = 1, BORN = 2, NONRES = 3 };
enum{ CHI2ONLY = 1, CHI2PENALTY = 2, ADAPTIVE = 3 };
enum{ NONE = 0, MLP1 = 1, MLP2 = 2, MLP3 = 3, CGLN = 4 };
enum{ SIG_0  =  0,
      SIG_S  =  1, SIG_T  =  2, SIG_P  =  3, SIG_E  =  4, SIG_F  =  5, SIG_G  =  6, SIG_H  =  7,
      SIG_CX =  8, SIG_CZ =  9, SIG_OX = 10, SIG_OZ = 11,
      ASY_S  = 12, ASY_T  = 13, ASY_P  = 14, ASY_E  = 15, ASY_F  = 16, ASY_G  = 17, ASY_H  = 18,
      ASY_CX = 19, ASY_CZ = 20, ASY_OX = 21, ASY_OZ = 22 };
static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;
static const Int_t LBINS  = 10;
static const Int_t FINITE = 256;
static const Int_t SOL = 16;
static const Int_t OBS = 23;
#endif
