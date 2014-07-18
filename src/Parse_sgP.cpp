#include "Parse_sgP.h"

//-----------------------------------------------------------------------------

void Parse_sgP()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaP, DsigmaP;
  FILE* File_sgP;

  printf("Loading sgP  data... ");
  File_sgP = fopen("data/sgP.txt", "r");

  sgP_bin = 0;
  while(!feof(File_sgP))
  {
    //Get beam energy
    if(fscanf(File_sgP, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgP(Energy);
    if(EnergyBin==sgP_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgP_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sgP, "%lf %lf %lf\n", &Theta, &sigmaP, &DsigmaP)==3)
    {
      sgP_val[EnergyBin][ThetaBin] = sigmaP;
      sgP_err[EnergyBin][ThetaBin] = DsigmaP;
      sgP_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaP!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgP);

    sgP_pts[EnergyBin] = ThetaBin;
    sgP_en[EnergyBin] = Energy;
    sgP_wt[EnergyBin] = Weight;
    sgP_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgP_bin)
      sgP_bin++;
  }

  fclose(File_sgP);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgP_bin; t++) n+=sgP_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgP_bin; t++) if(sgP_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgP_bin);
  for(Int_t e=0; e<sgP_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgP_en[e], sgP_pts[e]);
    for(Int_t th=0; th<sgP_pts[e]; th++)
      printf("%f %f %f\n", sgP_th[e][th], sgP_val[e][th], sgP_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgP()
{
  //Get energy bin for sigma0*P for given global energy
  Double_t Min = 1e38;
  Int_t eP = 0;

  for(Int_t e=0; e<sgP_bin; e++)
    if(fabs(sgP_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgP_en[e] - gEnergy);
      eP = e;
    }
  return eP;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgP(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eP = sgP_bin;

  for(Int_t e=0; e<sgP_bin; e++)
    if(sgP_en[e]==Energy)
      eP = e;

  return eP;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgP()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgP;
  Int_t eP = GetEnergyBin_sgP();

  //Calculate chi^2 for sigma0*P data
  ChiSq_sgP = 0.0;
  Omega = sgP_en[eP];
  for(Int_t th=0; th<sgP_pts[eP]; th++)
  {
    Theta = sgP_th[eP][th];
    Meas  = sgP_val[eP][th]*f_obs[SIG_P];
    Error = sgP_err[eP][th]*f_obs[SIG_P];
    Theo  = sigmaP(Theta, Omega);
    //printf("sgP: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgP+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgP_wt[eP]*ChiSq_sgP; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgP()
{
  Int_t eP = GetEnergyBin_sgP();
  Double_t Scale_sgP = (1.0*sgP_pts[eP])*(f_obs[SIG_P]-1.0)*(f_obs[SIG_P]-1.0)/(sgP_sy[eP]*sgP_sy[eP]);
  return Scale_sgP;
}

//-----------------------------------------------------------------------------
