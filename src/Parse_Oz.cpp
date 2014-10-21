#include "Parse_Oz.h"

//-----------------------------------------------------------------------------

void Parse_Oz()
{
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, Oz, DOz;
  FILE* File_Oz;

  printf("Loading   Oz data... ");
  File_Oz = fopen("data/Oz.txt", "r");

  Oz_bin = 0;
  while(!feof(File_Oz))
  {
    //Get beam energy
    if(fscanf(File_Oz, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_Oz(Energy);
    if(EnergyBin==Oz_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = Oz_pts[EnergyBin]; //..hence we append to the existing energy bin

    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Oz(File_Oz, &Theta, &Oz, &DOz)==3)
    {
      Oz_val[EnergyBin][ThetaBin] = Oz;
      Oz_err[EnergyBin][ThetaBin] = DOz;
      Oz_th[EnergyBin][ThetaBin]  = Theta;
      if(DOz!=0.0) ThetaBin++; //Accept only 'existing' data points
    }

    Oz_pts[EnergyBin] = ThetaBin;
    Oz_en[EnergyBin] = Energy;
    Oz_wt[EnergyBin] = Weight;
    Oz_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==Oz_bin)
      Oz_bin++;
  }

  fclose(File_Oz);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Oz_bin; t++) n+=Oz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Oz_bin; t++) if(Oz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Oz_bin);
  for(Int_t e=0; e<Oz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Oz_en[e], Oz_pts[e]);
    for(Int_t th=0; th<Oz_pts[e]; th++)
      printf("%f %f %f\n", Oz_th[e][th], Oz_val[e][th], Oz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_Oz()
{
  //Get energy bin for Oz for given global energy
  Double_t Min = 1e38;
  Int_t eO = 0;

  for(Int_t e=0; e<Oz_bin; e++)
    if(fabs(Oz_en[e] - gEnergy) < Min)
    {
      Min = fabs(Oz_en[e] - gEnergy);
      eO = e;
    }
  return eO;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_Oz(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eO = Oz_bin;

  for(Int_t e=0; e<Oz_bin; e++)
    if(Oz_en[e]==Energy)
      eO = e;

  return eO;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Oz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Oz;
  Int_t eO = GetEnergyBin_Oz();

  //Calculate chi^2 for Oz data
  ChiSq_Oz = 0.0;
  Omega = Oz_en[eO];
  for(Int_t th=0; th<Oz_pts[eO]; th++)
  {
    Theta = Oz_th[eO][th];
    Meas  = Oz_val[eO][th]*f_obs[ASY_OZ];
    Error = Oz_err[eO][th]*f_obs[ASY_OZ];
    Theo  = Oz(Theta, Omega);
    //printf("Oz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_Oz+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return Oz_wt[eO]*ChiSq_Oz; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_Oz()
{
  Int_t eO = GetEnergyBin_Oz();
  Double_t Scale_Oz = (1.0*Oz_pts[eO])*(f_obs[ASY_OZ]-1.0)*(f_obs[ASY_OZ]-1.0)/(Oz_sy[eO]*Oz_sy[eO]);
  return Scale_Oz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Oz(FILE* File_Oz, Double_t* Theta, Double_t* Oz, Double_t* DOz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Oz);
  return sscanf(Buffer, "%lf %lf %lf", Theta, Oz, DOz);
}

//-----------------------------------------------------------------------------
