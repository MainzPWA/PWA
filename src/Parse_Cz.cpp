#include "Parse_Cz.h"

//-----------------------------------------------------------------------------

void Parse_Cz()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, Cz, DCz;
  FILE* File_Cz;

  printf("Loading   Cz data... ");
  File_Cz = fopen("data/Cz.txt", "r");

  Cz_bin = 0;
  while(!feof(File_Cz))
  {
    //Get beam energy
    if(fscanf(File_Cz, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_Cz(Energy);
    if(EnergyBin==Cz_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = Cz_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_Cz, "%lf %lf %lf\n", &Theta, &Cz, &DCz)==3)
    {
      Cz_val[EnergyBin][ThetaBin] = Cz;
      Cz_err[EnergyBin][ThetaBin] = DCz;
      Cz_th[EnergyBin][ThetaBin]  = Theta;
      if((Cz!=0.0) && (DCz!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_Cz);

    Cz_pts[EnergyBin] = ThetaBin;
    Cz_en[EnergyBin] = Energy;
    Cz_wt[EnergyBin] = Weight;
    Cz_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==Cz_bin)
      Cz_bin++;
  }

  fclose(File_Cz);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Cz_bin; t++) n+=Cz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Cz_bin; t++) if(Cz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Cz_bin);
  for(Int_t e=0; e<Cz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Cz_en[e], Cz_pts[e]);
    for(Int_t th=0; th<Cz_pts[e]; th++)
      printf("%f %f %f\n", Cz_th[e][th], Cz_val[e][th], Cz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_Cz()
{
  //Get energy bin for Cz for given global energy
  Double_t Min = 1e38;
  Int_t eC = 0;

  for(Int_t e=0; e<Cz_bin; e++)
    if(fabs(Cz_en[e] - gEnergy) < Min)
    {
      Min = fabs(Cz_en[e] - gEnergy);
      eC = e;
    }
  return eC;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_Cz(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eC = Cz_bin;

  for(Int_t e=0; e<Cz_bin; e++)
    if(Cz_en[e]==Energy)
      eC = e;

  return eC;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Cz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Cz;
  Int_t eC = GetEnergyBin_Cz();

  //Calculate chi^2 for Cz data
  ChiSq_Cz = 0.0;
  Omega = Cz_en[eC];
  for(Int_t th=0; th<Cz_pts[eC]; th++)
  {
    Theta = Cz_th[eC][th];
    Meas  = Cz_val[eC][th]*f_obs[ASY_CZ];
    Error = Cz_err[eC][th]*f_obs[ASY_CZ];
    Theo  = Cz(Theta, Omega);
    //printf("Cz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_Cz+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return Cz_wt[eC]*ChiSq_Cz; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_Cz()
{
  Int_t eC = GetEnergyBin_Cz();
  Double_t Scale_Cz = (1.0*Cz_pts[eC])*(f_obs[ASY_CZ]-1.0)*(f_obs[ASY_CZ]-1.0)/(Cz_sy[eC]*Cz_sy[eC]);
  return Scale_Cz;
}

//-----------------------------------------------------------------------------
