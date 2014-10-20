#include "Parse_sgE.h"

//-----------------------------------------------------------------------------

void Parse_sgE()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaE, DsigmaE, Dummy;
  FILE* File_sgE;

  printf("Loading sgE  data... ");
  File_sgE = fopen("data/sgE.txt", "r");

  sgE_bin = 0;
  while(!feof(File_sgE))
  {
    //Get beam energy
    if(fscanf(File_sgE, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgE(Energy);
    if(EnergyBin==sgE_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgE_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sgE, "%lf %lf %lf\n", &Theta, &sigmaE, &DsigmaE, &Dummy)>=3)
    {
      sgE_val[EnergyBin][ThetaBin] = sigmaE;
      sgE_err[EnergyBin][ThetaBin] = DsigmaE;
      sgE_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaE!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgE);

    sgE_pts[EnergyBin] = ThetaBin;
    sgE_en[EnergyBin] = Energy;
    sgE_wt[EnergyBin] = Weight;
    sgE_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgE_bin)
      sgE_bin++;
  }

  fclose(File_sgE);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgE_bin; t++) n+=sgE_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgE_bin; t++) if(sgE_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgE_bin);
  for(Int_t e=0; e<sgE_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgE_en[e], sgE_pts[e]);
    for(Int_t th=0; th<sgE_pts[e]; th++)
      printf("%f %f %f\n", sgE_th[e][th], sgE_val[e][th], sgE_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgE()
{
  //Get energy bin for sigma0*E for given global energy
  Double_t Min = 1e38;
  Int_t eE = 0;

  for(Int_t e=0; e<sgE_bin; e++)
    if(fabs(sgE_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgE_en[e] - gEnergy);
      eE = e;
    }
  return eE;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgE(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eE = sgE_bin;

  for(Int_t e=0; e<sgE_bin; e++)
    if(sgE_en[e]==Energy)
      eE = e;

  return eE;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgE()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgE;
  Int_t eE = GetEnergyBin_sgE();

  //Calculate chi^2 for sigma0*E data
  ChiSq_sgE = 0.0;
  Omega = sgE_en[eE];
  for(Int_t th=0; th<sgE_pts[eE]; th++)
  {
    Theta = sgE_th[eE][th];
    Meas  = sgE_val[eE][th]*f_obs[SIG_E];
    Error = sgE_err[eE][th]*f_obs[SIG_E];
    Theo  = sigmaE(Theta, Omega);
    //printf("sgE: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgE+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgE_wt[eE]*ChiSq_sgE; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgE()
{
  Int_t eE = GetEnergyBin_sgE();
  Double_t Scale_sgE = (1.0*sgE_pts[eE])*(f_obs[SIG_E]-1.0)*(f_obs[SIG_E]-1.0)/(sgE_sy[eE]*sgE_sy[eE]);
  return Scale_sgE;
}

//-----------------------------------------------------------------------------
