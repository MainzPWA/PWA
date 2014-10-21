#include "Parse_sgG.h"

//-----------------------------------------------------------------------------

void Parse_sgG()
{
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaG, DsigmaG;
  FILE* File_sgG;

  printf("Loading sgG  data... ");
  File_sgG = fopen("data/sgG.txt", "r");

  sgG_bin = 0;
  while(!feof(File_sgG))
  {
    //Get beam energy
    if(fscanf(File_sgG, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgG(Energy);
    if(EnergyBin==sgG_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgG_pts[EnergyBin]; //..hence we append to the existing energy bin

    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgG(File_sgG, &Theta, &sigmaG, &DsigmaG)==3)
    {
      sgG_val[EnergyBin][ThetaBin] = sigmaG;
      sgG_err[EnergyBin][ThetaBin] = DsigmaG;
      sgG_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaG!=0.0) ThetaBin++; //Accept only 'existing' data points
    }

    sgG_pts[EnergyBin] = ThetaBin;
    sgG_en[EnergyBin] = Energy;
    sgG_wt[EnergyBin] = Weight;
    sgG_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgG_bin)
      sgG_bin++;
  }

  fclose(File_sgG);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgG_bin; t++) n+=sgG_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgG_bin; t++) if(sgG_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgG_bin);
  for(Int_t e=0; e<sgG_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgG_en[e], sgG_pts[e]);
    for(Int_t th=0; th<sgG_pts[e]; th++)
      printf("%f %f %f\n", sgG_th[e][th], sgG_val[e][th], sgG_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgG()
{
  //Get energy bin for sigma0*G for given global energy
  Double_t Min = 1e38;
  Int_t eG = 0;

  for(Int_t e=0; e<sgG_bin; e++)
    if(fabs(sgG_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgG_en[e] - gEnergy);
      eG = e;
    }
  return eG;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgG(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eG = sgG_bin;

  for(Int_t e=0; e<sgG_bin; e++)
    if(sgG_en[e]==Energy)
      eG = e;

  return eG;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgG()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgG;
  Int_t eG = GetEnergyBin_sgG();

  //Calculate chi^2 for sigma0*G data
  ChiSq_sgG = 0.0;
  Omega = sgG_en[eG];
  for(Int_t th=0; th<sgG_pts[eG]; th++)
  {
    Theta = sgG_th[eG][th];
    Meas  = sgG_val[eG][th]*f_obs[SIG_G];
    Error = sgG_err[eG][th]*f_obs[SIG_G];
    Theo  = sigmaG(Theta, Omega);
    //printf("sgG: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgG+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgG_wt[eG]*ChiSq_sgG; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgG()
{
  Int_t eG = GetEnergyBin_sgG();
  Double_t Scale_sgG = (1.0*sgG_pts[eG])*(f_obs[SIG_G]-1.0)*(f_obs[SIG_G]-1.0)/(sgG_sy[eG]*sgG_sy[eG]);
  return Scale_sgG;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgG(FILE* File_sgG, Double_t* Theta, Double_t* sigmaG, Double_t* DsigmaG)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgG);
  return sscanf(Buffer, "%lf %lf %lf", Theta, sigmaG, DsigmaG);
}

//-----------------------------------------------------------------------------
