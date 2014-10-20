#include "Parse_sgOx.h"

//-----------------------------------------------------------------------------

void Parse_sgOx()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaOx, DsigmaOx;
  FILE* File_sgOx;

  printf("Loading sgOx data... ");
  File_sgOx = fopen("data/sgOx.txt", "r");

  sgOx_bin = 0;
  while(!feof(File_sgOx))
  {
    //Get beam energy
    if(fscanf(File_sgOx, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgOx(Energy);
    if(EnergyBin==sgOx_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgOx_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sgOx, "%lf %lf %lf %*\n", &Theta, &sigmaOx, &DsigmaOx)==3)
    {
      sgOx_val[EnergyBin][ThetaBin] = sigmaOx;
      sgOx_err[EnergyBin][ThetaBin] = DsigmaOx;
      sgOx_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaOx!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgOx);

    sgOx_pts[EnergyBin] = ThetaBin;
    sgOx_en[EnergyBin] = Energy;
    sgOx_wt[EnergyBin] = Weight;
    sgOx_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgOx_bin)
      sgOx_bin++;
  }

  fclose(File_sgOx);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgOx_bin; t++) n+=sgOx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgOx_bin; t++) if(sgOx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgOx_bin);
  for(Int_t e=0; e<sgOx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgOx_en[e], sgOx_pts[e]);
    for(Int_t th=0; th<sgOx_pts[e]; th++)
      printf("%f %f %f\n", sgOx_th[e][th], sgOx_val[e][th], sgOx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgOx()
{
  //Get energy bin for sigmaOx for given global energy
  Double_t Min = 1e38;
  Int_t eO = 0;

  for(Int_t e=0; e<sgOx_bin; e++)
    if(fabs(sgOx_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgOx_en[e] - gEnergy);
      eO = e;
    }
  return eO;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgOx(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eO = sgOx_bin;

  for(Int_t e=0; e<sgOx_bin; e++)
    if(sgOx_en[e]==Energy)
      eO = e;

  return eO;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgOx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgOx;
  Int_t eO = GetEnergyBin_sgOx();

  //Calculate chi^2 for sigmaOx data
  ChiSq_sgOx = 0.0;
  Omega = sgOx_en[eO];
  for(Int_t th=0; th<sgOx_pts[eO]; th++)
  {
    Theta = sgOx_th[eO][th];
    Meas  = sgOx_val[eO][th]*f_obs[SIG_OX];
    Error = sgOx_err[eO][th]*f_obs[SIG_OX];
    Theo  = sigmaOx(Theta, Omega);
    //printf("sgOx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgOx+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgOx_wt[eO]*ChiSq_sgOx; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgOx()
{
  Int_t eO = GetEnergyBin_sgOx();
  Double_t Scale_sgOx = (1.0*sgOx_pts[eO])*(f_obs[SIG_OX]-1.0)*(f_obs[SIG_OX]-1.0)/(sgOx_sy[eO]*sgOx_sy[eO]);
  return Scale_sgOx;
}

//-----------------------------------------------------------------------------
