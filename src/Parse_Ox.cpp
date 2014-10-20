#include "Parse_Ox.h"

//-----------------------------------------------------------------------------

void Parse_Ox()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, Ox, DOx, Dummy;
  FILE* File_Ox;

  printf("Loading   Ox data... ");
  File_Ox = fopen("data/Ox.txt", "r");

  Ox_bin = 0;
  while(!feof(File_Ox))
  {
    //Get beam energy
    if(fscanf(File_Ox, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_Ox(Energy);
    if(EnergyBin==Ox_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = Ox_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_Ox, "%lf %lf %lf\n", &Theta, &Ox, &DOx, &Dummy)>=3)
    {
      Ox_val[EnergyBin][ThetaBin] = Ox;
      Ox_err[EnergyBin][ThetaBin] = DOx;
      Ox_th[EnergyBin][ThetaBin]  = Theta;
      if(DOx!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_Ox);

    Ox_pts[EnergyBin] = ThetaBin;
    Ox_en[EnergyBin] = Energy;
    Ox_wt[EnergyBin] = Weight;
    Ox_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==Ox_bin)
      Ox_bin++;
  }

  fclose(File_Ox);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Ox_bin; t++) n+=Ox_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Ox_bin; t++) if(Ox_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Ox_bin);
  for(Int_t e=0; e<Ox_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Ox_en[e], Ox_pts[e]);
    for(Int_t th=0; th<Ox_pts[e]; th++)
      printf("%f %f %f\n", Ox_th[e][th], Ox_val[e][th], Ox_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_Ox()
{
  //Get energy bin for Ox for given global energy
  Double_t Min = 1e38;
  Int_t eO = 0;

  for(Int_t e=0; e<Ox_bin; e++)
    if(fabs(Ox_en[e] - gEnergy) < Min)
    {
      Min = fabs(Ox_en[e] - gEnergy);
      eO = e;
    }
  return eO;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_Ox(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eO = Ox_bin;

  for(Int_t e=0; e<Ox_bin; e++)
    if(Ox_en[e]==Energy)
      eO = e;

  return eO;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Ox()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Ox;
  Int_t eO = GetEnergyBin_Ox();

  //Calculate chi^2 for Ox data
  ChiSq_Ox = 0.0;
  Omega = Ox_en[eO];
  for(Int_t th=0; th<Ox_pts[eO]; th++)
  {
    Theta = Ox_th[eO][th];
    Meas  = Ox_val[eO][th]*f_obs[ASY_OX];
    Error = Ox_err[eO][th]*f_obs[ASY_OX];
    Theo  = Ox(Theta, Omega);
    //printf("Ox: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_Ox+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return Ox_wt[eO]*ChiSq_Ox; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_Ox()
{
  Int_t eO = GetEnergyBin_Ox();
  Double_t Scale_Ox = (1.0*Ox_pts[eO])*(f_obs[ASY_OX]-1.0)*(f_obs[ASY_OX]-1.0)/(Ox_sy[eO]*Ox_sy[eO]);
  return Scale_Ox;
}

//-----------------------------------------------------------------------------
