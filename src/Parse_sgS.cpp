#include "Parse_sgS.h"

//-----------------------------------------------------------------------------

void Parse_sgS()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaS, DsigmaS;
  FILE* File_sgS;

  printf("Loading sgS  data... ");
  File_sgS = fopen("data/sgS.txt", "r");

  sgS_bin = 0;
  while(!feof(File_sgS))
  {
    //Get beam energy
    if(fscanf(File_sgS, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgS(Energy);
    if(EnergyBin==sgS_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgS_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sgS, "%lf %lf %lf\n", &Theta, &sigmaS, &DsigmaS)==3)
    {
      sgS_val[EnergyBin][ThetaBin] = sigmaS;
      sgS_err[EnergyBin][ThetaBin] = DsigmaS;
      sgS_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaS!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgS);

    sgS_pts[EnergyBin] = ThetaBin;
    sgS_en[EnergyBin] = Energy;
    sgS_wt[EnergyBin] = Weight;
    sgS_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgS_bin)
      sgS_bin++;
  }

  fclose(File_sgS);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgS_bin; t++) n+=sgS_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgS_bin; t++) if(sgS_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgS_bin);
  for(Int_t e=0; e<sgS_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgS_en[e], sgS_pts[e]);
    for(Int_t th=0; th<sgS_pts[e]; th++)
      printf("%f %f %f\n", sgS_th[e][th], sgS_val[e][th], sgS_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgS()
{
  //Get energy bin for sigma0*Sigma for given global energy
  Double_t Min = 1e38;
  Int_t eS = 0;

  for(Int_t e=0; e<sgS_bin; e++)
    if(fabs(sgS_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgS_en[e] - gEnergy);
      eS = e;
    }
  return eS;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgS(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eS = sgS_bin;

  for(Int_t e=0; e<sgS_bin; e++)
    if(sgS_en[e]==Energy)
      eS = e;

  return eS;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgS()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgS;
  Int_t eS = GetEnergyBin_sgS();

  //Calculate chi^2 for sigma0*Sigma data
  ChiSq_sgS = 0.0;
  Omega = sgS_en[eS];
  for(Int_t th=0; th<sgS_pts[eS]; th++)
  {
    Theta = sgS_th[eS][th];
    Meas  = sgS_val[eS][th]*f_obs[SIG_S];
    Error = sgS_err[eS][th]*f_obs[SIG_S];
    Theo  = sigmaS(Theta, Omega);
    //printf("sgS: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgS+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgS_wt[eS]*ChiSq_sgS; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgS()
{
  Int_t eS = GetEnergyBin_sgS();
  Double_t Scale_sgS = (1.0*sgS_pts[eS])*(f_obs[SIG_S]-1.0)*(f_obs[SIG_S]-1.0)/(sgS_sy[eS]*sgS_sy[eS]);
  return Scale_sgS;
}

//-----------------------------------------------------------------------------
