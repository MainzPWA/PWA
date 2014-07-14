#include "Parse_H.h"

//-----------------------------------------------------------------------------

void Parse_H()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, H, DH;
  FILE* File_H;

  printf("Loading   H data... ");
  File_H = fopen("data/H.txt", "r");

  H_bin = 0;
  while(!feof(File_H))
  {
    //Get beam energy
    if(fscanf(File_H, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_H(Energy);
    if(EnergyBin==H_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = H_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_H, "%lf %lf %lf\n", &Theta, &H, &DH)==3)
    {
      H_val[EnergyBin][ThetaBin] = H;
      H_err[EnergyBin][ThetaBin] = DH;
      H_th[EnergyBin][ThetaBin]  = Theta;
      if((H!=0.0) && (DH!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_H);

    H_pts[EnergyBin] = ThetaBin;
    H_en[EnergyBin] = Energy;
    H_wt[EnergyBin] = Weight;
    H_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==H_bin)
      H_bin++;
  }

  fclose(File_H);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<H_bin; t++) n+=H_pts[t];
  Int_t m = 0; for(Int_t t=0; t<H_bin; t++) if(H_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", H_bin);
  for(Int_t e=0; e<H_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, H_en[e], H_pts[e]);
    for(Int_t th=0; th<H_pts[e]; th++)
      printf("%f %f %f\n", H_th[e][th], H_val[e][th], H_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_H()
{
  //Get energy bin for H for given global energy
  Double_t Min = 1e38;
  Int_t eH = 0;

  for(Int_t e=0; e<H_bin; e++)
    if(fabs(H_en[e] - gEnergy) < Min)
    {
      Min = fabs(H_en[e] - gEnergy);
      eH = e;
    }
  return eH;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_H(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eH = H_bin;

  for(Int_t e=0; e<H_bin; e++)
    if(H_en[e]==Energy)
      eH = e;

  return eH;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_H()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_H;
  Int_t eH = GetEnergyBin_H();

  //Calculate chi^2 for H data
  ChiSq_H = 0.0;
  Omega = H_en[eH];
  for(Int_t th=0; th<H_pts[eH]; th++)
  {
    Theta = H_th[eH][th];
    Meas  = H_val[eH][th]*f_obs[ASY_H];
    Error = H_err[eH][th]*f_obs[ASY_H];
    Theo  = H(Theta, Omega);
    //printf("H: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_H+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return H_wt[eH]*ChiSq_H; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_H()
{
  Int_t eH = GetEnergyBin_H();
  Double_t Scale_H = (1.0*H_pts[eH])*(f_obs[ASY_H]-1.0)*(f_obs[ASY_H]-1.0)/(H_sy[eH]*H_sy[eH]);
  return Scale_H;
}

//-----------------------------------------------------------------------------
