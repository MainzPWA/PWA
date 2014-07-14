#include "Parse_S.h"

//-----------------------------------------------------------------------------

void Parse_S()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, Sigma, DSigma;
  FILE* File_S;

  printf("Loading   S data... ");
  File_S = fopen("data/S.txt", "r");

  S_bin = 0;
  while(!feof(File_S))
  {
    //Get beam energy
    if(fscanf(File_S, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_S(Energy);
    if(EnergyBin==S_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = S_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_S, "%lf %lf %lf\n", &Theta, &Sigma, &DSigma)==3)
    {
      S_val[EnergyBin][ThetaBin] = Sigma;
      S_err[EnergyBin][ThetaBin] = DSigma;
      S_th[EnergyBin][ThetaBin]  = Theta;
      if((Sigma!=0.0) && (DSigma!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_S);

    S_pts[EnergyBin] = ThetaBin;
    S_en[EnergyBin] = Energy;
    S_wt[EnergyBin] = Weight;
    S_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==S_bin)
      S_bin++;
  }

  fclose(File_S);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<S_bin; t++) n+=S_pts[t];
  Int_t m = 0; for(Int_t t=0; t<S_bin; t++) if(S_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", S_bin);
  for(Int_t e=0; e<S_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, S_en[e], S_pts[e]);
    for(Int_t th=0; th<S_pts[e]; th++)
      printf("%f %f %f\n", S_th[e][th], S_val[e][th], S_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_S()
{
  //Get energy bin for Sigma for given global energy
  Double_t Min = 1e38;
  Int_t eS = 0;

  for(Int_t e=0; e<S_bin; e++)
    if(fabs(S_en[e] - gEnergy) < Min)
    {
      Min = fabs(S_en[e] - gEnergy);
      eS = e;
    }
  return eS;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_S(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eS = S_bin;

  for(Int_t e=0; e<S_bin; e++)
    if(S_en[e]==Energy)
      eS = e;

  return eS;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_S()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_S;
  Int_t eS = GetEnergyBin_S();

  //Calculate chi^2 for Sigma data
  ChiSq_S = 0.0;
  Omega = S_en[eS];
  for(Int_t th=0; th<S_pts[eS]; th++)
  {
    Theta = S_th[eS][th];
    Meas  = S_val[eS][th]*f_obs[ASY_S];
    Error = S_err[eS][th]*f_obs[ASY_S];
    Theo  = S(Theta, Omega);
    //printf("S: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_S+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return S_wt[eS]*ChiSq_S; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_S()
{
  Int_t eS = GetEnergyBin_S();
  Double_t Scale_S = (1.0*S_pts[eS])*(f_obs[ASY_S]-1.0)*(f_obs[ASY_S]-1.0)/(S_sy[eS]*S_sy[eS]);
  return Scale_S;
}

//-----------------------------------------------------------------------------
