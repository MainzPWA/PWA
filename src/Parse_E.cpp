#include "Parse_E.h"

//-----------------------------------------------------------------------------

void Parse_E()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, E, DE;
  FILE* File_E;

  printf("Loading   E  data... ");
  File_E = fopen("data/E.txt", "r");

  E_bin = 0;
  while(!feof(File_E))
  {
    //Get beam energy
    if(fscanf(File_E, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_E(Energy);
    if(EnergyBin==E_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = E_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_E, "%lf %lf %lf\n", &Theta, &E, &DE)==3)
    {
      E_val[EnergyBin][ThetaBin] = E;
      E_err[EnergyBin][ThetaBin] = DE;
      E_th[EnergyBin][ThetaBin]  = Theta;
      if(DE!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_E);

    E_pts[EnergyBin] = ThetaBin;
    E_en[EnergyBin] = Energy;
    E_wt[EnergyBin] = Weight;
    E_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==E_bin)
      E_bin++;
  }

  fclose(File_E);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<E_bin; t++) n+=E_pts[t];
  Int_t m = 0; for(Int_t t=0; t<E_bin; t++) if(E_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", E_bin);
  for(Int_t e=0; e<E_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, E_en[e], E_pts[e]);
    for(Int_t th=0; th<E_pts[e]; th++)
      printf("%f %f %f\n", E_th[e][th], E_val[e][th], E_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_E()
{
  //Get energy bin for E for given global energy
  Double_t Min = 1e38;
  Int_t eE = 0;

  for(Int_t e=0; e<E_bin; e++)
    if(fabs(E_en[e] - gEnergy) < Min)
    {
      Min = fabs(E_en[e] - gEnergy);
      eE = e;
    }
  return eE;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_E(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eE = E_bin;

  for(Int_t e=0; e<E_bin; e++)
    if(E_en[e]==Energy)
      eE = e;

  return eE;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_E()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_E;
  Int_t eE = GetEnergyBin_E();

  //Calculate chi^2 for E data
  ChiSq_E = 0.0;
  Omega = E_en[eE];
  for(Int_t th=0; th<E_pts[eE]; th++)
  {
    Theta = E_th[eE][th];
    Meas  = E_val[eE][th]*f_obs[ASY_E];
    Error = E_err[eE][th]*f_obs[ASY_E];
    Theo  = E(Theta, Omega);
    //printf("E: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_E+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return E_wt[eE]*ChiSq_E; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_E()
{
  Int_t eE = GetEnergyBin_E();
  Double_t Scale_E = (1.0*E_pts[eE])*(f_obs[ASY_E]-1.0)*(f_obs[ASY_E]-1.0)/(E_sy[eE]*E_sy[eE]);
  return Scale_E;
}

//-----------------------------------------------------------------------------
