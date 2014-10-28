#include "Parse_sgCx.h"

//-----------------------------------------------------------------------------

void Parse_sgCx()
{
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaCx, DsigmaCx;
  FILE* File_sgCx;

  printf("Loading sgCx data... ");
  File_sgCx = fopen("data/sgCx.txt", "r");

  sgCx_bin = 0;
  while(!feof(File_sgCx))
  {
    //Get beam energy
    if(fscanf(File_sgCx, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgCx(Energy);
    if(EnergyBin==sgCx_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgCx_pts[EnergyBin]; //..hence we append to the existing energy bin

    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgCx(File_sgCx, &Theta, &sigmaCx, &DsigmaCx)==3)
    {
      sgCx_val[EnergyBin][ThetaBin] = sigmaCx;
      sgCx_err[EnergyBin][ThetaBin] = DsigmaCx;
      sgCx_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaCx!=0.0) ThetaBin++; //Accept only 'existing' data points
    }

    sgCx_pts[EnergyBin] = ThetaBin;
    sgCx_en[EnergyBin] = Energy;
    sgCx_wt[EnergyBin] = Weight;
    sgCx_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgCx_bin)
      sgCx_bin++;
  }

  fclose(File_sgCx);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgCx_bin; t++) n+=sgCx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgCx_bin; t++) if(sgCx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgCx_bin);
  for(Int_t e=0; e<sgCx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgCx_en[e], sgCx_pts[e]);
    for(Int_t th=0; th<sgCx_pts[e]; th++)
      printf("%f %f %f\n", sgCx_th[e][th], sgCx_val[e][th], sgCx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgCx()
{
  //Get energy bin for sigmaCx for given global energy
  Double_t Min = 1e38;
  Int_t eC = 0;

  for(Int_t e=0; e<sgCx_bin; e++)
    if(fabs(sgCx_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgCx_en[e] - gEnergy);
      eC = e;
    }
  return eC;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgCx(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eC = sgCx_bin;

  for(Int_t e=0; e<sgCx_bin; e++)
    if(sgCx_en[e]==Energy)
      eC = e;

  return eC;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgCx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgCx;
  Int_t eC = GetEnergyBin_sgCx();

  //Calculate chi^2 for sigmaCx data
  ChiSq_sgCx = 0.0;
  Omega = sgCx_en[eC];
  for(Int_t th=0; th<sgCx_pts[eC]; th++)
  {
    Theta = sgCx_th[eC][th];
    Meas  = sgCx_val[eC][th]*f_obs[SIG_CX];
    Error = sgCx_err[eC][th]*f_obs[SIG_CX];
    Theo  = sigmaCx(Theta, Omega);
    //printf("sgCx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgCx+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgCx_wt[eC]*ChiSq_sgCx; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgCx()
{
  Int_t eC = GetEnergyBin_sgCx();
  Double_t Scale_sgCx = (1.0*sgCx_pts[eC])*(f_obs[SIG_CX]-1.0)*(f_obs[SIG_CX]-1.0)/(sgCx_sy[eC]*sgCx_sy[eC]);
  return Scale_sgCx;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgCx(FILE* File_sgCx, Double_t* Theta, Double_t* sigmaCx, Double_t* DsigmaCx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgCx);
  return sscanf(Buffer, "%lf %lf %lf", Theta, sigmaCx, DsigmaCx);
}

//-----------------------------------------------------------------------------
