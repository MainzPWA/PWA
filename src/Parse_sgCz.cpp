#include "Parse_sgCz.h"

//-----------------------------------------------------------------------------

void Parse_sgCz()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaCz, DsigmaCz;
  FILE* File_sgCz;

  printf("Loading sgCz data... ");
  File_sgCz = fopen("data/sgCz.txt", "r");

  sgCz_bin = 0;
  while(!feof(File_sgCz))
  {
    //Get beam energy
    if(fscanf(File_sgCz, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgCz(Energy);
    if(EnergyBin==sgCz_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgCz_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sgCz, "%lf %lf %lf %*\n", &Theta, &sigmaCz, &DsigmaCz)==3)
    {
      sgCz_val[EnergyBin][ThetaBin] = sigmaCz;
      sgCz_err[EnergyBin][ThetaBin] = DsigmaCz;
      sgCz_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaCz!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgCz);

    sgCz_pts[EnergyBin] = ThetaBin;
    sgCz_en[EnergyBin] = Energy;
    sgCz_wt[EnergyBin] = Weight;
    sgCz_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgCz_bin)
      sgCz_bin++;
  }

  fclose(File_sgCz);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgCz_bin; t++) n+=sgCz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgCz_bin; t++) if(sgCz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgCz_bin);
  for(Int_t e=0; e<sgCz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgCz_en[e], sgCz_pts[e]);
    for(Int_t th=0; th<sgCz_pts[e]; th++)
      printf("%f %f %f\n", sgCz_th[e][th], sgCz_val[e][th], sgCz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgCz()
{
  //Get energy bin for sigmaCz for given global energy
  Double_t Min = 1e38;
  Int_t eC = 0;

  for(Int_t e=0; e<sgCz_bin; e++)
    if(fabs(sgCz_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgCz_en[e] - gEnergy);
      eC = e;
    }
  return eC;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgCz(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eC = sgCz_bin;

  for(Int_t e=0; e<sgCz_bin; e++)
    if(sgCz_en[e]==Energy)
      eC = e;

  return eC;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgCz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgCz;
  Int_t eC = GetEnergyBin_sgCz();

  //Calculate chi^2 for sigmaCz data
  ChiSq_sgCz = 0.0;
  Omega = sgCz_en[eC];
  for(Int_t th=0; th<sgCz_pts[eC]; th++)
  {
    Theta = sgCz_th[eC][th];
    Meas  = sgCz_val[eC][th]*f_obs[SIG_CZ];
    Error = sgCz_err[eC][th]*f_obs[SIG_CZ];
    Theo  = sigmaCz(Theta, Omega);
    //printf("sgCz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgCz+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgCz_wt[eC]*ChiSq_sgCz; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgCz()
{
  Int_t eC = GetEnergyBin_sgCz();
  Double_t Scale_sgCz = (1.0*sgCz_pts[eC])*(f_obs[SIG_CZ]-1.0)*(f_obs[SIG_CZ]-1.0)/(sgCz_sy[eC]*sgCz_sy[eC]);
  return Scale_sgCz;
}

//-----------------------------------------------------------------------------
