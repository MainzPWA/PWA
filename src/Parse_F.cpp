#include "Parse_F.h"

//-----------------------------------------------------------------------------

void Parse_F()
{
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, F, DF;
  FILE* File_F;

  printf("Loading   F  data... ");
  File_F = fopen("data/F.txt", "r");

  F_bin = 0;
  while(!feof(File_F))
  {
    //Get beam energy
    if(fscanf(File_F, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_F(Energy);
    if(EnergyBin==F_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = F_pts[EnergyBin]; //..hence we append to the existing energy bin

    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_F(File_F, &Theta, &F, &DF)==3)
    {
      F_val[EnergyBin][ThetaBin] = F;
      F_err[EnergyBin][ThetaBin] = DF;
      F_th[EnergyBin][ThetaBin]  = Theta;
      if(DF!=0.0) ThetaBin++; //Accept only 'existing' data points
    }

    F_pts[EnergyBin] = ThetaBin;
    F_en[EnergyBin] = Energy;
    F_wt[EnergyBin] = Weight;
    F_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==F_bin)
      F_bin++;
  }

  fclose(File_F);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<F_bin; t++) n+=F_pts[t];
  Int_t m = 0; for(Int_t t=0; t<F_bin; t++) if(F_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", F_bin);
  for(Int_t e=0; e<F_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, F_en[e], F_pts[e]);
    for(Int_t th=0; th<F_pts[e]; th++)
      printf("%f %f %f\n", F_th[e][th], F_val[e][th], F_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_F()
{
  //Get energy bin for F for given global energy
  Double_t Min = 1e38;
  Int_t eF = 0;

  for(Int_t e=0; e<F_bin; e++)
    if(fabs(F_en[e] - gEnergy) < Min)
    {
      Min = fabs(F_en[e] - gEnergy);
      eF = e;
    }
  return eF;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_F(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eF = F_bin;

  for(Int_t e=0; e<F_bin; e++)
    if(F_en[e]==Energy)
      eF = e;

  return eF;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_F()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_F;
  Int_t eF = GetEnergyBin_F();

  //Calculate chi^2 for F data
  ChiSq_F = 0.0;
  Omega = F_en[eF];
  for(Int_t th=0; th<F_pts[eF]; th++)
  {
    Theta = F_th[eF][th];
    Meas  = F_val[eF][th]*f_obs[ASY_F];
    Error = F_err[eF][th]*f_obs[ASY_F];
    Theo  = F(Theta, Omega);
    //printf("F: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_F+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return F_wt[eF]*ChiSq_F; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_F()
{
  Int_t eF = GetEnergyBin_F();
  Double_t Scale_F = (1.0*F_pts[eF])*(f_obs[ASY_F]-1.0)*(f_obs[ASY_F]-1.0)/(F_sy[eF]*F_sy[eF]);
  return Scale_F;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_F(FILE* File_F, Double_t* Theta, Double_t* F, Double_t* DF)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_F);
  return sscanf(Buffer, "%lf %lf %lf", Theta, F, DF);
}

//-----------------------------------------------------------------------------
