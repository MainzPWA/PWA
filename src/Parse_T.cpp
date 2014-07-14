#include "Parse_T.h"

//-----------------------------------------------------------------------------

void Parse_T()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, T, DT;
  FILE* File_T;

  printf("Loading   T data... ");
  File_T = fopen("data/T.txt", "r");

  T_bin = 0;
  while(!feof(File_T))
  {
    //Get beam energy
    if(fscanf(File_T, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_T(Energy);
    if(EnergyBin==T_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = T_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_T, "%lf %lf %lf\n", &Theta, &T, &DT)==3)
    {
      T_val[EnergyBin][ThetaBin] = T;
      T_err[EnergyBin][ThetaBin] = DT;
      T_th[EnergyBin][ThetaBin]  = Theta;
      if((T!=0.0) && (DT!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_T);

    T_pts[EnergyBin] = ThetaBin;
    T_en[EnergyBin] = Energy;
    T_wt[EnergyBin] = Weight;
    T_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==T_bin)
      T_bin++;
  }

  fclose(File_T);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<T_bin; t++) n+=T_pts[t];
  Int_t m = 0; for(Int_t t=0; t<T_bin; t++) if(T_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", T_bin);
  for(Int_t e=0; e<T_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, T_en[e], T_pts[e]);
    for(Int_t th=0; th<T_pts[e]; th++)
      printf("%f %f %f\n", T_th[e][th], T_val[e][th], T_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_T()
{
  //Get energy bin for Sigma for given global energy
  Double_t Min = 1e38;
  Int_t eT = 0;

  for(Int_t e=0; e<T_bin; e++)
    if(fabs(T_en[e] - gEnergy) < Min)
    {
      Min = fabs(T_en[e] - gEnergy);
      eT = e;
    }
  return eT;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_T(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eT = T_bin;

  for(Int_t e=0; e<T_bin; e++)
    if(T_en[e]==Energy)
      eT = e;

  return eT;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_T()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_T;
  Int_t eT = GetEnergyBin_T();

  //Calculate chi^2 for T data
  ChiSq_T = 0.0;
  Omega = T_en[eT];
  for(Int_t th=0; th<T_pts[eT]; th++)
  {
    Theta = T_th[eT][th];
    Meas  = T_val[eT][th]*f_obs[ASY_T];
    Error = T_err[eT][th]*f_obs[ASY_T];
    Theo  = T(Theta, Omega);
    //printf("T: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_T+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return T_wt[eT]*ChiSq_T; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_T()
{
  Int_t eT = GetEnergyBin_T();
  Double_t Scale_T = (1.0*T_pts[eT])*(f_obs[ASY_T]-1.0)*(f_obs[ASY_T]-1.0)/(T_sy[eT]*T_sy[eT]);
  return Scale_T;
}

//-----------------------------------------------------------------------------
