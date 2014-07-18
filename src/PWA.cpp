#include "Constants.h"
#include "Globals.h"
#include "Multipoles.h"
#include "Fitter.h"
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
#include "PWA.h"
#include "Version.h"
#include "Build.h"

//-----------------------------------------------------------------------------

void PlotFilename(Char_t* Buffer, Int_t p, Int_t s)
{
  if(p==0) sprintf(Buffer, "plots.%d/ReE0p.root", s);
  if(p==1) sprintf(Buffer, "plots.%d/ImE0p.root", s);
  if(p==2) sprintf(Buffer, "plots.%d/ReE1p.root", s);
  if(p==3) sprintf(Buffer, "plots.%d/ImE1p.root", s);
  if(p==4) sprintf(Buffer, "plots.%d/ReM1p.root", s);
  if(p==5) sprintf(Buffer, "plots.%d/ImM1p.root", s);
  if(p==6) sprintf(Buffer, "plots.%d/ReM1m.root", s);
  if(p==7) sprintf(Buffer, "plots.%d/ImM1m.root", s);

  if(p > 7)
  {
    if(p%8==0) sprintf(Buffer, "plots.%d/ReE%dp.root", s, p/8 + 1);
    if(p%8==1) sprintf(Buffer, "plots.%d/ImE%dp.root", s, p/8 + 1);
    if(p%8==2) sprintf(Buffer, "plots.%d/ReM%dp.root", s, p/8 + 1);
    if(p%8==3) sprintf(Buffer, "plots.%d/ImM%dp.root", s, p/8 + 1);
    if(p%8==4) sprintf(Buffer, "plots.%d/ReE%dm.root", s, p/8 + 1);
    if(p%8==5) sprintf(Buffer, "plots.%d/ImE%dm.root", s, p/8 + 1);
    if(p%8==6) sprintf(Buffer, "plots.%d/ReM%dm.root", s, p/8 + 1);
    if(p%8==7) sprintf(Buffer, "plots.%d/ImM%dm.root", s, p/8 + 1);
  }
}

//-----------------------------------------------------------------------------

void PlotTitle(Char_t* Buffer, Int_t p)
{
  if(p==0) sprintf(Buffer, "ReE0p");
  if(p==1) sprintf(Buffer, "ImE0p");
  if(p==2) sprintf(Buffer, "ReE1p");
  if(p==3) sprintf(Buffer, "ImE1p");
  if(p==4) sprintf(Buffer, "ReM1p");
  if(p==5) sprintf(Buffer, "ImM1p");
  if(p==6) sprintf(Buffer, "ReM1m");
  if(p==7) sprintf(Buffer, "ImM1m");

  if(p > 7)
  {
    if(p%8==0) sprintf(Buffer, "ReE%dp", p/8 + 1);
    if(p%8==1) sprintf(Buffer, "ImE%dp", p/8 + 1);
    if(p%8==2) sprintf(Buffer, "ReM%dp", p/8 + 1);
    if(p%8==3) sprintf(Buffer, "ImM%dp", p/8 + 1);
    if(p%8==4) sprintf(Buffer, "ReE%dm", p/8 + 1);
    if(p%8==5) sprintf(Buffer, "ImE%dm", p/8 + 1);
    if(p%8==6) sprintf(Buffer, "ReM%dm", p/8 + 1);
    if(p%8==7) sprintf(Buffer, "ImM%dm", p/8 + 1);
  }
}

//-----------------------------------------------------------------------------

void Plot(Int_t s)
{
  TCanvas* Plot;
  TGraphErrors* Fit;
  TGraph* Chi;
  TGraph* Pen;
  TGraph* Model;
  Char_t Name[256];
  Char_t Title[256];
  Char_t Filename[256];

  //Plot s,p,...,L_max wave multipoles
  for(Int_t t=0; t<L_MAX*8; t++)
  {
    PlotTitle(Title, t);
    PlotFilename(Filename, t, s);
    Plot = new TCanvas(Title, Title);

    sprintf(Name, "%s_fit", Title);
    Fit = new TGraphErrors(Fit_pts[s], Fit_en[s], Fit_val[s][t], NULL, Fit_err[s][t]);
    Fit->SetTitle(Title);
    Fit->SetName(Name);
    Fit->SetMarkerSize(0.8);
    if(t%2)
    {
      Fit->SetMarkerStyle(24);
      Fit->SetMarkerColor(kBlue+2);
      Fit->SetLineColor(kBlue+2);
    }
    else
    {
      Fit->SetMarkerStyle(20);
      Fit->SetMarkerColor(kRed+2);
      Fit->SetLineColor(kRed+2);
    }
    Fit->Draw("APZ");
    Fit->GetXaxis()->SetTitle("#omega / MeV");

    sprintf(Name, "%s_model", Title);
    Model = new TGraph(Model_pts, Model_en, Model_val[t]);
    Model->SetTitle(Title);
    Model->SetName(Name);
    if(t%2)
      Model->SetLineColor(kBlue);
    else
      Model->SetLineColor(kRed);
    Model->SetLineWidth(2);
    Model->Draw("L");

    gDirectory->SaveObjectAs(Plot, Filename, "q");
    delete Plot;
  }

  //Plot chi^2
  Plot = new TCanvas("Chi2", "Chi2");
  Chi = new TGraph(Fit_pts[s], Fit_en[s], Fit_chi[s]);
  Chi->SetTitle("Chi2");
  Chi->SetName("Chi2");
  Chi->SetMarkerSize(0.8);
  Chi->SetMarkerStyle(20);
  Chi->Draw("APZ");
  Chi->GetXaxis()->SetTitle("#omega / MeV");
  sprintf(Filename, "plots.%1d/Chi2.root", s);
  gDirectory->SaveObjectAs(Plot, Filename, "q");
  delete Plot;

  //Plot penalty
  Plot = new TCanvas("Penalty", "Penalty");
  Pen = new TGraph(Fit_pts[s], Fit_en[s], Fit_pen[s]);
  Pen->SetTitle("Penalty");
  Pen->SetName("Penalty");
  Pen->SetMarkerSize(0.8);
  Pen->SetMarkerStyle(20);
  Pen->Draw("APZ");
  Pen->GetXaxis()->SetTitle("#omega / MeV");
  sprintf(Filename, "plots.%1d/Penalty.root", s);
  gDirectory->SaveObjectAs(Plot, Filename, "q");
  delete Plot;
}

//-----------------------------------------------------------------------------

void Print()
{
  //Select MAID energy closest to cross section energy
  Int_t eM = GetEnergyBin_maid();
  Double_t DeviationSq;

  printf("\n------------------------------------------------------------------------------------\n");
  printf("omega = %7.3f MeV, W = %8.3f MeV\n", gEnergy, W_cm(gEnergy));
  printf("Mlp  |  Fit (Chi^2/NDF = %6.3f, NDF = %3d)       |  Model               |  Chi'^2\n", ChiSq()/NDF(), NDF());
  printf("-----+--------------------------------------------+----------------------+----------\n");

  if(ONLY_CROSS_S || ONLY_CROSS_F) //Handle special cases for s,p wave imaginary parts in threshold fits
  {
    //This is rather dirty, as it works destructive on the loaded model parameters.
    //However, the overwritten values are never used anyway...
    maid_Ep[0][eM] = TComplex(maid_Ep[0][eM].Re(), ImE0p());
    maid_Ep[1][eM] = TComplex(maid_Ep[1][eM].Re(), 0.0);
    maid_Mp[1][eM] = TComplex(maid_Mp[1][eM].Re(), 0.0);
    maid_Mm[1][eM] = TComplex(maid_Mm[1][eM].Re(), 0.0);
  }
  DeviationSq = (Ep[0]-maid_Ep[0][eM]).Rho2()/DEp[0].Rho2();
  printf("E0+  |  (%7.3f +- %5.3f) + (%7.3f +- %5.3f)i  |  %7.3f + %7.3fi  |  %8.4f\n",
          Ep[0].Re(), DEp[0].Re(), Ep[0].Im(), DEp[0].Im(), maid_Ep[0][eM].Re(), maid_Ep[0][eM].Im(), DeviationSq);
  DeviationSq = (Ep[1]-maid_Ep[1][eM]).Rho2()/DEp[1].Rho2();
  printf("E1+  |  (%7.3f +- %5.3f) + (%7.3f +- %5.3f)i  |  %7.3f + %7.3fi  |  %8.4f\n",
          Ep[1].Re(), DEp[1].Re(), Ep[1].Im(), DEp[1].Im(), maid_Ep[1][eM].Re(), maid_Ep[1][eM].Im(), DeviationSq);
  DeviationSq = (Mp[1]-maid_Mp[1][eM]).Rho2()/DMp[1].Rho2();
  printf("M1+  |  (%7.3f +- %5.3f) + (%7.3f +- %5.3f)i  |  %7.3f + %7.3fi  |  %8.4f\n",
          Mp[1].Re(), DMp[1].Re(), Mp[1].Im(), DMp[1].Im(), maid_Mp[1][eM].Re(), maid_Mp[1][eM].Im(), DeviationSq);
  DeviationSq = (Mm[1]-maid_Mm[1][eM]).Rho2()/DMm[1].Rho2();
  printf("M1-  |  (%7.3f +- %5.3f) + (%7.3f +- %5.3f)i  |  %7.3f + %7.3fi  |  %8.4f\n",
          Mm[1].Re(), DMm[1].Re(), Mm[1].Im(), DMm[1].Im(), maid_Mm[1][eM].Re(), maid_Mm[1][eM].Im(), DeviationSq);

  for(Int_t l=2; l<L_MAX+1; l++)
  {
    if(FIX_EP[l]) DeviationSq = 0.0; else DeviationSq = (Ep[l]-maid_Ep[l][eM]).Rho2()/DEp[l].Rho2();
    printf("E%d+  |  (%7.3f +- %5.3f) + (%7.3f +- %5.3f)i  |  %7.3f + %7.3fi  |  %8.4f\n",
            l, Ep[l].Re(), DEp[l].Re(), Ep[l].Im(), DEp[l].Im(), maid_Ep[l][eM].Re(), maid_Ep[l][eM].Im(), DeviationSq);
    if(FIX_EM[l]) DeviationSq = 0.0; else DeviationSq = (Em[l]-maid_Em[l][eM]).Rho2()/DEm[l].Rho2();
    printf("E%d-  |  (%7.3f +- %5.3f) + (%7.3f +- %5.3f)i  |  %7.3f + %7.3fi  |  %8.4f\n",
            l, Em[l].Re(), DEm[l].Re(), Em[l].Im(), DEm[l].Im(), maid_Em[l][eM].Re(), maid_Em[l][eM].Im(), DeviationSq);
    if(FIX_MP[l]) DeviationSq = 0.0; else DeviationSq = (Mp[l]-maid_Mp[l][eM]).Rho2()/DMp[l].Rho2();
    printf("M%d+  |  (%7.3f +- %5.3f) + (%7.3f +- %5.3f)i  |  %7.3f + %7.3fi  |  %8.4f\n",
            l, Mp[l].Re(), DMp[l].Re(), Mp[l].Im(), DMp[l].Im(), maid_Mp[l][eM].Re(), maid_Mp[l][eM].Im(), DeviationSq);
    if(FIX_MM[l]) DeviationSq = 0.0; else DeviationSq = (Mm[l]-maid_Mm[l][eM]).Rho2()/DMm[l].Rho2();
    printf("M%d-  |  (%7.3f +- %5.3f) + (%7.3f +- %5.3f)i  |  %7.3f + %7.3fi  |  %8.4f\n",
            l, Mm[l].Re(), DMm[l].Re(), Mm[l].Im(), DMm[l].Im(), maid_Mm[l][eM].Re(), maid_Mm[l][eM].Im(), DeviationSq);
  }

  if(FIX_SCALES)
    printf("-----+--------------------------------------------+----------------------+----------\n");
  else
  {
    printf("-----+------------------+------+------------------+----------------------+----------\n");
    printf("sgA  |  Scale           |  A   |  Scale           |\n");
    printf("-----+------------------+------+------------------+---------------------------------\n");
    printf("sg0  |  %5.3f +- %5.3f  |      |                  |\n", f_obs[SIG_0], Df_obs[SIG_0]);
    printf("sgS  |  %5.3f +- %5.3f  | ",  f_obs[SIG_S],  Df_obs[SIG_S]);
    printf(" S   |  %5.3f +- %5.3f  |\n", f_obs[ASY_S],  Df_obs[ASY_S]);
    printf("sgT  |  %5.3f +- %5.3f  | ",  f_obs[SIG_T],  Df_obs[SIG_T]);
    printf(" T   |  %5.3f +- %5.3f  |\n", f_obs[ASY_T],  Df_obs[ASY_T]);
    printf("sgP  |  %5.3f +- %5.3f  | ",  f_obs[SIG_P],  Df_obs[SIG_P]);
    printf(" P   |  %5.3f +- %5.3f  |\n", f_obs[ASY_P],  Df_obs[ASY_P]);
    printf("sgE  |  %5.3f +- %5.3f  | ",  f_obs[SIG_E],  Df_obs[SIG_E]);
    printf(" E   |  %5.3f +- %5.3f  |\n", f_obs[ASY_E],  Df_obs[ASY_E]);
    printf("sgF  |  %5.3f +- %5.3f  | ",  f_obs[SIG_F],  Df_obs[SIG_F]);
    printf(" F   |  %5.3f +- %5.3f  |\n", f_obs[ASY_F],  Df_obs[ASY_F]);
    printf("sgG  |  %5.3f +- %5.3f  | ",  f_obs[SIG_G],  Df_obs[SIG_G]);
    printf(" G   |  %5.3f +- %5.3f  |\n", f_obs[ASY_G],  Df_obs[ASY_G]);
    printf("sgH  |  %5.3f +- %5.3f  | ",  f_obs[SIG_H],  Df_obs[SIG_H]);
    printf(" H   |  %5.3f +- %5.3f  |\n", f_obs[ASY_H],  Df_obs[ASY_H]);
    printf("sgCx |  %5.3f +- %5.3f  | ",  f_obs[SIG_CX], Df_obs[SIG_CX]);
    printf(" Cx  |  %5.3f +- %5.3f  |\n", f_obs[ASY_CX], Df_obs[ASY_CX]);
    printf("sgCz |  %5.3f +- %5.3f  | ",  f_obs[SIG_CZ], Df_obs[SIG_CZ]);
    printf(" Cz  |  %5.3f +- %5.3f  |\n", f_obs[ASY_CZ], Df_obs[ASY_CZ]);
    printf("sgOx |  %5.3f +- %5.3f  | ",  f_obs[SIG_OX], Df_obs[SIG_OX]);
    printf(" Ox  |  %5.3f +- %5.3f  |\n", f_obs[ASY_OX], Df_obs[ASY_OX]);
    printf("sgOz |  %5.3f +- %5.3f  | ",  f_obs[SIG_OZ], Df_obs[SIG_OZ]);
    printf(" Oz  |  %5.3f +- %5.3f  |\n", f_obs[ASY_OZ], Df_obs[ASY_OZ]);
    printf("------------------------+------+------------------+---------------------------------\n");
  }

  if(PRINT_PENALTY)
  {
    printf("Chi^2 = %7.3f, Penalty = %7.3f, Scale = %7.3f, NPts = %2d, NPar = %2d, NDF = %2d\n",
           ChiSq(), Penalty(), Scale(), NPts(), NPar(), NDF());
    printf("------------------------------------------------------------------------------------\n");
  }
}

//-----------------------------------------------------------------------------

void PWA()
{
  for(Int_t s=0; s<SOLUTIONS; s++) Fit_pts[s] = 0;
  Model_pts = 0;

  for(Int_t e0=0; e0<sg0_bin; e0++)
  {
    //Set global energy to cross section energy
    gEnergy = sg0_en[e0];
    //Fit only within selected energy range
    if((gEnergy < MIN_ENERGY) || (gEnergy > MAX_ENERGY)) continue;
    //Fit only if enough data available
    if(NDF() < 1) continue;
    //Fit only if cross section data available
    if(sg0_pts[e0]==0) continue;

    //Perform fit
    Fit();

    //Print out fit results for best solution...
    Print();
    //...and store fit result plots for all solutions to .root file
    for(Int_t s=0; s<SOLUTIONS; s++) Plot(s);
  }
  printf("\n");
}

//-----------------------------------------------------------------------------

void Init()
{
  FILE* Config;
  Char_t Buffer[1024];
  Double_t Double;
  Bool_t Bool;
  Int_t Int;

  //Configuration options for fitting process
  FIX_EP_PHASE[0] = false;
  FIX_EP_PHASE[1] = false;
  FIX_MP_PHASE[1] = false;
  FIX_MM_PHASE[1] = false;
  for(Int_t l=2; l<LBINS; l++)
  {
    FIX_EP[l] = true; FIX_EP_PHASE[l] = false;
    FIX_EM[l] = true; FIX_EM_PHASE[l] = false;
    FIX_MP[l] = true; FIX_MP_PHASE[l] = false;
    FIX_MM[l] = true; FIX_MM_PHASE[l] = false;
  }
  FIX_SCALES    = false;
  FIX_IM_E0P    = false;
  ONLY_CROSS_S  = false;
  ONLY_CROSS_F  = false;
  PRINT_PENALTY = false;
  ERROR_MODE    = ADAPTIVE;
  PENALTY_MODE  = MLP1;
  D_WAVES       = MODEL;
  L_MAX         = 2;
  PENALTY[MLP1] = 1.0;
  PENALTY[MLP2] = 1.0;
  SCALING       = 0.1;
  MIN_ENERGY    = 144.7;
  MAX_ENERGY    = 420.0;
  VARIATION     = 0.20;
  BETA          = 3.43;
  ITERATIONS    = 16;
  SOLUTIONS     = 1;
  MASS_MESON    = MASS_PIZERO;
  MASS_INITIAL  = MASS_PROTON;
  MASS_FINAL    = MASS_PROTON;
  MASS2_MESON   = MASS2_PIZERO;
  MASS2_INITIAL = MASS2_PROTON;
  MASS2_FINAL   = MASS2_PROTON;
  THRESHOLD     = W_thres();

  Config = fopen("PWA.cfg", "r");
  if(!Config) return;

  while(!feof(Config))
  {
    fgets(Buffer, sizeof(Buffer), Config);
    if(sscanf(Buffer, "ITERATIONS %d", &Int)==1) ITERATIONS = Int;
    if(sscanf(Buffer, "SOLUTIONS %d", &Int)==1) SOLUTIONS = Int;
    if(sscanf(Buffer, "ERROR_MODE %d", &Int)==1) ERROR_MODE = Int;
    if(sscanf(Buffer, "PENALTY_MODE %d", &Int)==1) PENALTY_MODE = Int;
    if(sscanf(Buffer, "L_MAX %d", &Int)==1) L_MAX = Int;
    if(sscanf(Buffer, "D_WAVES %d", &Int)==1) D_WAVES = Int;
    if(sscanf(Buffer, "FIX_E%dP %d", &Int, &Bool)==2) FIX_EP[Int] = Bool;
    if(sscanf(Buffer, "FIX_M%dP %d", &Int, &Bool)==2) FIX_MP[Int] = Bool;
    if(sscanf(Buffer, "FIX_E%dM %d", &Int, &Bool)==2) FIX_EM[Int] = Bool;
    if(sscanf(Buffer, "FIX_M%dM %d", &Int, &Bool)==2) FIX_MM[Int] = Bool;
    if(sscanf(Buffer, "FIX_IM_E0P %d", &Bool)==1) FIX_IM_E0P = Bool;
    if(sscanf(Buffer, "FIX_E%dP_PHASE %d", &Int, &Bool)==2) FIX_EP_PHASE[Int] = Bool;
    if(sscanf(Buffer, "FIX_E%dM_PHASE %d", &Int, &Bool)==2) FIX_EM_PHASE[Int] = Bool;
    if(sscanf(Buffer, "FIX_M%dP_PHASE %d", &Int, &Bool)==2) FIX_MP_PHASE[Int] = Bool;
    if(sscanf(Buffer, "FIX_M%dM_PHASE %d", &Int, &Bool)==2) FIX_MM_PHASE[Int] = Bool;
    if(sscanf(Buffer, "FIX_SCALES %d", &Bool)==1) FIX_SCALES = Bool;
    if(sscanf(Buffer, "ONLY_CROSS_S %d", &Bool)==1) ONLY_CROSS_S = Bool;
    if(sscanf(Buffer, "ONLY_CROSS_F %d", &Bool)==1) ONLY_CROSS_F = Bool;
    if(sscanf(Buffer, "PRINT_PENALTY %d", &Bool)==1) PRINT_PENALTY = Bool;
    if(sscanf(Buffer, "PENALTY_MLP1 %lf", &Double)==1) PENALTY[MLP1] = Double;
    if(sscanf(Buffer, "PENALTY_MLP2 %lf", &Double)==1) PENALTY[MLP2] = Double;
    if(sscanf(Buffer, "SCALING %lf", &Double)==1) SCALING = Double;
    if(sscanf(Buffer, "BETA %lf", &Double)==1) BETA = Double;
    if(sscanf(Buffer, "VARIATION %lf", &Double)==1) VARIATION = Double;
    if(sscanf(Buffer, "MIN_ENERGY %lf", &Double)==1) MIN_ENERGY = Double;
    if(sscanf(Buffer, "MAX_ENERGY %lf", &Double)==1) MAX_ENERGY = Double;
    if(sscanf(Buffer, "MASS_MESON %lf", &Double)==1) MASS_MESON = Double;
    if(sscanf(Buffer, "MASS_INITIAL %lf", &Double)==1) MASS_INITIAL = Double;
    if(sscanf(Buffer, "MASS_FINAL %lf", &Double)==1) MASS_FINAL = Double;
  }
  fclose(Config);
  MASS2_MESON   = MASS_MESON*MASS_MESON;
  MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;
  MASS2_FINAL   = MASS_FINAL*MASS_FINAL;
  THRESHOLD     = W_thres();

  //Check some nonsenical parameter combinations
  if(MAX_ENERGY < MIN_ENERGY) //Invalid energy range
  {
    printf("Error: Invalid energy range.\n");
    exit(0);
  }
  if(MIN_ENERGY < THRESHOLD) //Sub-threshold energy
  {
    printf("Error: MIN_ENERGY below threshold.\n");
    exit(0);
  }
  if(L_MAX+1 > LBINS) //Angular momentum too big
  {
    printf("Error: L_MAX must be less than %d.\n", LBINS);
    exit(0);
  }
  if(SOLUTIONS > SOL) //Too many solutions requested
  {
    printf("Error: SOLUTIONS must be less than %d.\n", SOL+1);
    exit(0);
  }
  if(SOLUTIONS > ITERATIONS) //Too many solutions requested
  {
    printf("Error: SOLUTIONS must not be larger than ITERATIONS.\n");
    exit(0);
  }
  if(ONLY_CROSS_S && ONLY_CROSS_F) //Either fit sigma0,S or sigma0,F but not both
  {
    printf("Error: Illegal combination of ONLY_CROSS_S and ONLY_CROSS_F.\n");
    exit(0);
  }
  if(FIX_IM_E0P && !(ONLY_CROSS_S || ONLY_CROSS_F)) //Fixed imaginary part of E0+ is only allowed for limited observable set
  {
    printf("Error: FIX_IM_E0P only allowed for ONLY_CROSS_S or ONLY_CROSS_F.\n");
    exit(0);
  }
  if((ONLY_CROSS_S || ONLY_CROSS_F) && (FIX_EP_PHASE[0] || FIX_EP_PHASE[1] || FIX_MP_PHASE[1] || FIX_MM_PHASE[1]))
  {
    printf("Error: ONLY_CROSS_S or ONLY_CROSS_F must not be used with explicit phase fixing.\n");
    exit(0);
  }
  Bool = FIX_EP_PHASE[0] || FIX_EP_PHASE[1] || FIX_MP_PHASE[1] || FIX_MM_PHASE[1];
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    Bool = Bool || (FIX_EP[l] || FIX_EM[l] || FIX_MP[l] || FIX_MM[l]); //Check for fixed multipoles
    Bool = Bool || (FIX_EP_PHASE[l] || FIX_EM_PHASE[l] || FIX_MP_PHASE[l] || FIX_MM_PHASE[l]); //Check for fixed phases
  }
  if(!Bool && !(ONLY_CROSS_S || ONLY_CROSS_F)) //In general, a phase should be fixed. However, with the penalty fitting, it might work even will all phases open
    printf("Warning: Fixing at least one phase might be necessary.\n");

  printf("------------------------------------------------------------------------------------\n");
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    printf("FIX_E%dP %1d\n", l, FIX_EP[l]);
    printf("FIX_M%dP %1d\n", l, FIX_MP[l]);
    printf("FIX_E%dM %1d\n", l, FIX_EM[l]);
    printf("FIX_M%dM %1d\n", l, FIX_MM[l]);
  }
  printf("FIX_IM_E0P %1d\n", FIX_IM_E0P);
  printf("FIX_E0P_PHASE %1d\n", FIX_EP_PHASE[0]);
  printf("FIX_E1P_PHASE %1d\n", FIX_EP_PHASE[1]);
  printf("FIX_M1P_PHASE %1d\n", FIX_MP_PHASE[1]);
  printf("FIX_M1M_PHASE %1d\n", FIX_MM_PHASE[1]);
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    printf("FIX_E%dP_PHASE %1d\n", l, FIX_EP_PHASE[l]);
    printf("FIX_M%dP_PHASE %1d\n", l, FIX_MP_PHASE[l]);
    printf("FIX_E%dM_PHASE %1d\n", l, FIX_EM_PHASE[l]);
    printf("FIX_M%dM_PHASE %1d\n", l, FIX_MM_PHASE[l]);
  }
  printf("FIX_SCALES %1d\n", FIX_SCALES);
  printf("L_MAX %1d\n", L_MAX);
  printf("MIN_ENERGY %6.1f\n", MIN_ENERGY);
  printf("MAX_ENERGY %6.1f\n", MAX_ENERGY);
  printf("D_WAVES %1d\n", D_WAVES);
  printf("ONLY_CROSS_S %1d\n", ONLY_CROSS_S);
  printf("ONLY_CROSS_F %1d\n", ONLY_CROSS_F);
  printf("PENALTY_MLP1 %6.3f\n", PENALTY[MLP1]);
  printf("PENALTY_MLP2 %6.3f\n", PENALTY[MLP2]);
  printf("SCALING %6.3f\n", SCALING);
  printf("ERROR_MODE %1d\n", ERROR_MODE);
  printf("PENALTY_MODE %1d\n", PENALTY_MODE);
  printf("PRINT_PENALTY %1d\n", PRINT_PENALTY);
  printf("ITERATIONS %4d\n", ITERATIONS);
  printf("SOLUTIONS %2d\n", SOLUTIONS);
  printf("BETA %4.2f\n", BETA);
  printf("VARIATION %5.3f\n", VARIATION);
  printf("MASS_MESON %5.3f\n", MASS_MESON);
  printf("MASS_INITIAL %5.3f\n", MASS_INITIAL);
  printf("MASS_FINAL %5.3f\n", MASS_FINAL);
  printf("------------------------------------------------------------------------------------\n");
}

//-----------------------------------------------------------------------------

void Load()
{
  Parse_sg0();
  Parse_sgS();
  Parse_sgT();
  Parse_sgP();
  Parse_sgE();
  Parse_sgF();
  Parse_sgG();
  Parse_sgH();
  Parse_sgCx();
  Parse_sgCz();
  Parse_sgOx();
  Parse_sgOz();
  Parse_S();
  Parse_T();
  Parse_P();
  Parse_E();
  Parse_F();
  Parse_G();
  Parse_H();
  Parse_Cx();
  Parse_Cz();
  Parse_Ox();
  Parse_Oz();
  printf("------------------------------------------------------------------------------------\n");
  Parse_MAID();
  printf("------------------------------------------------------------------------------------\n");

  return;

  //Debug output
  for(Int_t e0=0; e0<sg0_bin; e0++)
    printf("%f: %d\n", sg0_en[e0], sg0_pts[e0]);
}

//-----------------------------------------------------------------------------

int main(int argc, char **argv)
{
  printf("\nPWA multipole fitter, version %.1f, build %d\n", VERSION, BUILD);

  gRandom->SetSeed(0);
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gMinuit = NULL;
  setbuf(stdout, NULL);

  Init();
  Load();
  PWA();
}

//-----------------------------------------------------------------------------

