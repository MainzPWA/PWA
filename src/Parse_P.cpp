#include "Parse_P.h"

//-----------------------------------------------------------------------------

void Parse_P()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, P, DP;
  FILE* File_P;

  printf("Loading   P  data... ");
  File_P = fopen("data/P.txt", "r");

  P_bin = 0;
  while(!feof(File_P))
  {
    //Get beam energy
    if(fscanf(File_P, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_P(Energy);
    if(EnergyBin==P_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = P_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_P, "%lf %lf %lf\n", &Theta, &P, &DP)==3)
    {
      P_val[EnergyBin][ThetaBin] = P;
      P_err[EnergyBin][ThetaBin] = DP;
      P_th[EnergyBin][ThetaBin]  = Theta;
      if(DP!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_P);

    P_pts[EnergyBin] = ThetaBin;
    P_en[EnergyBin] = Energy;
    P_wt[EnergyBin] = Weight;
    P_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==P_bin)
      P_bin++;
  }

  fclose(File_P);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<P_bin; t++) n+=P_pts[t];
  Int_t m = 0; for(Int_t t=0; t<P_bin; t++) if(P_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", P_bin);
  for(Int_t e=0; e<P_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, P_en[e], P_pts[e]);
    for(Int_t th=0; th<P_pts[e]; th++)
      printf("%f %f %f\n", P_th[e][th], P_val[e][th], P_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_P()
{
  //Get energy bin for P for given global energy
  Double_t Min = 1e38;
  Int_t eP = 0;

  for(Int_t e=0; e<P_bin; e++)
    if(fabs(P_en[e] - gEnergy) < Min)
    {
      Min = fabs(P_en[e] - gEnergy);
      eP = e;
    }
  return eP;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_P(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eP = P_bin;

  for(Int_t e=0; e<P_bin; e++)
    if(P_en[e]==Energy)
      eP = e;

  return eP;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_P()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_P;
  Int_t eP = GetEnergyBin_P();

  //Calculate chi^2 for P data
  ChiSq_P = 0.0;
  Omega = P_en[eP];
  for(Int_t th=0; th<P_pts[eP]; th++)
  {
    Theta = P_th[eP][th];
    Meas  = P_val[eP][th]*f_obs[ASY_P];
    Error = P_err[eP][th]*f_obs[ASY_P];
    Theo  = P(Theta, Omega);
    //printf("P: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_P+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return P_wt[eP]*ChiSq_P; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_P()
{
  Int_t eP = GetEnergyBin_P();
  Double_t Scale_P = (1.0*P_pts[eP])*(f_obs[ASY_P]-1.0)*(f_obs[ASY_P]-1.0)/(P_sy[eP]*P_sy[eP]);
  return Scale_P;
}

//-----------------------------------------------------------------------------
