#include "Parse_sgF.h"

//-----------------------------------------------------------------------------

void Parse_sgF()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaF, DsigmaF;
  FILE* File_sgF;

  printf("Loading sgF data... ");
  File_sgF = fopen("data/sgF.txt", "r");

  sgF_bin = 0;
  while(!feof(File_sgF))
  {
    //Get beam energy
    if(fscanf(File_sgF, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgF(Energy);
    if(EnergyBin==sgF_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgF_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sgF, "%lf %lf %lf\n", &Theta, &sigmaF, &DsigmaF)==3)
    {
      sgF_val[EnergyBin][ThetaBin] = sigmaF;
      sgF_err[EnergyBin][ThetaBin] = DsigmaF;
      sgF_th[EnergyBin][ThetaBin]  = Theta;
      if((sigmaF!=0.0) && (DsigmaF!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgF);

    sgF_pts[EnergyBin] = ThetaBin;
    sgF_en[EnergyBin] = Energy;
    sgF_wt[EnergyBin] = Weight;
    sgF_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgF_bin)
      sgF_bin++;
  }

  fclose(File_sgF);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgF_bin; t++) n+=sgF_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgF_bin; t++) if(sgF_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgF_bin);
  for(Int_t e=0; e<sgF_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgF_en[e], sgF_pts[e]);
    for(Int_t th=0; th<sgF_pts[e]; th++)
      printf("%f %f %f\n", sgF_th[e][th], sgF_val[e][th], sgF_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgF()
{
  //Get energy bin for sigma0*F for given global energy
  Double_t Min = 1e38;
  Int_t eF = 0;

  for(Int_t e=0; e<sgF_bin; e++)
    if(fabs(sgF_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgF_en[e] - gEnergy);
      eF = e;
    }
  return eF;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgF(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eF = sgF_bin;

  for(Int_t e=0; e<sgF_bin; e++)
    if(sgF_en[e]==Energy)
      eF = e;

  return eF;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgF()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgF;
  Int_t eF = GetEnergyBin_sgF();

  //Calculate chi^2 for sigma0*F data
  ChiSq_sgF = 0.0;
  Omega = sgF_en[eF];
  for(Int_t th=0; th<sgF_pts[eF]; th++)
  {
    Theta = sgF_th[eF][th];
    Meas  = sgF_val[eF][th]*f_obs[SIG_F];
    Error = sgF_err[eF][th]*f_obs[SIG_F];
    Theo  = sigmaF(Theta, Omega);
    //printf("sgF: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgF+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgF_wt[eF]*ChiSq_sgF; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgF()
{
  Int_t eF = GetEnergyBin_sgF();
  Double_t Scale_sgF = (1.0*sgF_pts[eF])*(f_obs[SIG_F]-1.0)*(f_obs[SIG_F]-1.0)/(sgF_sy[eF]*sgF_sy[eF]);
  return Scale_sgF;
}

//-----------------------------------------------------------------------------
