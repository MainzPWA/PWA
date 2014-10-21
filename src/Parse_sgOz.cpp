#include "Parse_sgOz.h"

//-----------------------------------------------------------------------------

void Parse_sgOz()
{
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaOz, DsigmaOz;
  FILE* File_sgOz;

  printf("Loading sgOz data... ");
  File_sgOz = fopen("data/sgOz.txt", "r");

  sgOz_bin = 0;
  while(!feof(File_sgOz))
  {
    //Get beam energy
    if(fscanf(File_sgOz, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgOz(Energy);
    if(EnergyBin==sgOz_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgOz_pts[EnergyBin]; //..hence we append to the existing energy bin

    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgOz(File_sgOz, &Theta, &sigmaOz, &DsigmaOz)==3)
    {
      sgOz_val[EnergyBin][ThetaBin] = sigmaOz;
      sgOz_err[EnergyBin][ThetaBin] = DsigmaOz;
      sgOz_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaOz!=0.0) ThetaBin++; //Accept only 'existing' data points
    }

    sgOz_pts[EnergyBin] = ThetaBin;
    sgOz_en[EnergyBin] = Energy;
    sgOz_wt[EnergyBin] = Weight;
    sgOz_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgOz_bin)
      sgOz_bin++;
  }

  fclose(File_sgOz);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgOz_bin; t++) n+=sgOz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgOz_bin; t++) if(sgOz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgOz_bin);
  for(Int_t e=0; e<sgOz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgOz_en[e], sgOz_pts[e]);
    for(Int_t th=0; th<sgOz_pts[e]; th++)
      printf("%f %f %f\n", sgOz_th[e][th], sgOz_val[e][th], sgOz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgOz()
{
  //Get energy bin for sigmaOz for given global energy
  Double_t Min = 1e38;
  Int_t eO = 0;

  for(Int_t e=0; e<sgOz_bin; e++)
    if(fabs(sgOz_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgOz_en[e] - gEnergy);
      eO = e;
    }
  return eO;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgOz(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eO = sgOz_bin;

  for(Int_t e=0; e<sgOz_bin; e++)
    if(sgOz_en[e]==Energy)
      eO = e;

  return eO;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgOz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgOz;
  Int_t eO = GetEnergyBin_sgOz();

  //Calculate chi^2 for sigmaOz data
  ChiSq_sgOz = 0.0;
  Omega = sgOz_en[eO];
  for(Int_t th=0; th<sgOz_pts[eO]; th++)
  {
    Theta = sgOz_th[eO][th];
    Meas  = sgOz_val[eO][th]*f_obs[SIG_OZ];
    Error = sgOz_err[eO][th]*f_obs[SIG_OZ];
    Theo  = sigmaOz(Theta, Omega);
    //printf("sgOz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgOz+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgOz_wt[eO]*ChiSq_sgOz; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgOz()
{
  Int_t eO = GetEnergyBin_sgOz();
  Double_t Scale_sgOz = (1.0*sgOz_pts[eO])*(f_obs[SIG_OZ]-1.0)*(f_obs[SIG_OZ]-1.0)/(sgOz_sy[eO]*sgOz_sy[eO]);
  return Scale_sgOz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgOz(FILE* File_sgOz, Double_t* Theta, Double_t* sigmaOz, Double_t* DsigmaOz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgOz);
  return sscanf(Buffer, "%lf %lf %lf", Theta, sigmaOz, DsigmaOz);
}

//-----------------------------------------------------------------------------
