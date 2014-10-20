#include "Parse_Cx.h"

//-----------------------------------------------------------------------------

void Parse_Cx()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, Cx, DCx, Dummy;
  FILE* File_Cx;

  printf("Loading   Cx data... ");
  File_Cx = fopen("data/Cx.txt", "r");

  Cx_bin = 0;
  while(!feof(File_Cx))
  {
    //Get beam energy
    if(fscanf(File_Cx, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_Cx(Energy);
    if(EnergyBin==Cx_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = Cx_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_Cx, "%lf %lf %lf\n", &Theta, &Cx, &DCx, &Dummy)>=3)
    {
      Cx_val[EnergyBin][ThetaBin] = Cx;
      Cx_err[EnergyBin][ThetaBin] = DCx;
      Cx_th[EnergyBin][ThetaBin]  = Theta;
      if(DCx!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_Cx);

    Cx_pts[EnergyBin] = ThetaBin;
    Cx_en[EnergyBin] = Energy;
    Cx_wt[EnergyBin] = Weight;
    Cx_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==Cx_bin)
      Cx_bin++;
  }

  fclose(File_Cx);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Cx_bin; t++) n+=Cx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Cx_bin; t++) if(Cx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Cx_bin);
  for(Int_t e=0; e<Cx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Cx_en[e], Cx_pts[e]);
    for(Int_t th=0; th<Cx_pts[e]; th++)
      printf("%f %f %f\n", Cx_th[e][th], Cx_val[e][th], Cx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_Cx()
{
  //Get energy bin for Cx for given global energy
  Double_t Min = 1e38;
  Int_t eC = 0;

  for(Int_t e=0; e<Cx_bin; e++)
    if(fabs(Cx_en[e] - gEnergy) < Min)
    {
      Min = fabs(Cx_en[e] - gEnergy);
      eC = e;
    }
  return eC;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_Cx(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eC = Cx_bin;

  for(Int_t e=0; e<Cx_bin; e++)
    if(Cx_en[e]==Energy)
      eC = e;

  return eC;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Cx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Cx;
  Int_t eC = GetEnergyBin_Cx();

  //Calculate chi^2 for Cx data
  ChiSq_Cx = 0.0;
  Omega = Cx_en[eC];
  for(Int_t th=0; th<Cx_pts[eC]; th++)
  {
    Theta = Cx_th[eC][th];
    Meas  = Cx_val[eC][th]*f_obs[ASY_CX];
    Error = Cx_err[eC][th]*f_obs[ASY_CX];
    Theo  = Cx(Theta, Omega);
    //printf("Cx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_Cx+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return Cx_wt[eC]*ChiSq_Cx; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_Cx()
{
  Int_t eC = GetEnergyBin_Cx();
  Double_t Scale_Cx = (1.0*Cx_pts[eC])*(f_obs[ASY_CX]-1.0)*(f_obs[ASY_CX]-1.0)/(Cx_sy[eC]*Cx_sy[eC]);
  return Scale_Cx;
}

//-----------------------------------------------------------------------------
