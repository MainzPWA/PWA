#include "Parse_sgH.h"

//-----------------------------------------------------------------------------

void Parse_sgH()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaH, DsigmaH;
  FILE* File_sgH;

  printf("Loading sgH  data... ");
  File_sgH = fopen("data/sgH.txt", "r");

  sgH_bin = 0;
  while(!feof(File_sgH))
  {
    //Get beam energy
    if(fscanf(File_sgH, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgH(Energy);
    if(EnergyBin==sgH_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgH_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sgH, "%lf %lf %lf\n", &Theta, &sigmaH, &DsigmaH)==3)
    {
      sgH_val[EnergyBin][ThetaBin] = sigmaH;
      sgH_err[EnergyBin][ThetaBin] = DsigmaH;
      sgH_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaH!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgH);

    sgH_pts[EnergyBin] = ThetaBin;
    sgH_en[EnergyBin] = Energy;
    sgH_wt[EnergyBin] = Weight;
    sgH_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgH_bin)
      sgH_bin++;
  }

  fclose(File_sgH);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgH_bin; t++) n+=sgH_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgH_bin; t++) if(sgH_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgH_bin);
  for(Int_t e=0; e<sgH_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgH_en[e], sgH_pts[e]);
    for(Int_t th=0; th<sgH_pts[e]; th++)
      printf("%f %f %f\n", sgH_th[e][th], sgH_val[e][th], sgH_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgH()
{
  //Get energy bin for sigma0*H for given global energy
  Double_t Min = 1e38;
  Int_t eH = 0;

  for(Int_t e=0; e<sgH_bin; e++)
    if(fabs(sgH_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgH_en[e] - gEnergy);
      eH = e;
    }
  return eH;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgH(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eH = sgH_bin;

  for(Int_t e=0; e<sgH_bin; e++)
    if(sgH_en[e]==Energy)
      eH = e;

  return eH;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgH()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgH;
  Int_t eH = GetEnergyBin_sgH();

  //Calculate chi^2 for sigma0*H data
  ChiSq_sgH = 0.0;
  Omega = sgH_en[eH];
  for(Int_t th=0; th<sgH_pts[eH]; th++)
  {
    Theta = sgH_th[eH][th];
    Meas  = sgH_val[eH][th]*f_obs[SIG_H];
    Error = sgH_err[eH][th]*f_obs[SIG_H];
    Theo  = sigmaH(Theta, Omega);
    //printf("sgH: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgH+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgH_wt[eH]*ChiSq_sgH; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgH()
{
  Int_t eH = GetEnergyBin_sgH();
  Double_t Scale_sgH = (1.0*sgH_pts[eH])*(f_obs[SIG_H]-1.0)*(f_obs[SIG_H]-1.0)/(sgH_sy[eH]*sgH_sy[eH]);
  return Scale_sgH;
}

//-----------------------------------------------------------------------------
