#include "Parse_sg0.h"

//-----------------------------------------------------------------------------

void Parse_sg0()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigma0, Dsigma0;
  FILE* File_sg0;

  printf("Loading sg0  data... ");
  File_sg0 = fopen("data/sg0.txt", "r");

  sg0_bin = 0;
  while(!feof(File_sg0))
  {
    //Get beam energy
    if(fscanf(File_sg0, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sg0(Energy);
    if(EnergyBin==sg0_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sg0_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sg0, "%lf %lf %lf\n", &Theta, &sigma0, &Dsigma0)==3)
    {
      sg0_val[EnergyBin][ThetaBin] = sigma0;
      sg0_err[EnergyBin][ThetaBin] = Dsigma0;
      sg0_th[EnergyBin][ThetaBin]  = Theta;
      if(Dsigma0!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sg0);

    sg0_pts[EnergyBin] = ThetaBin;
    sg0_en[EnergyBin] = Energy;
    sg0_wt[EnergyBin] = Weight;
    sg0_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sg0_bin)
      sg0_bin++;
  }

  fclose(File_sg0);
  Sort_sg0(0, sg0_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sg0_bin; t++) n+=sg0_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sg0_bin; t++) if(sg0_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sg0_bin);
  for(Int_t e=0; e<sg0_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sg0_en[e], sg0_pts[e]);
    for(Int_t th=0; th<sg0_pts[e]; th++)
      printf("%f %f %f\n", sg0_th[e][th], sg0_val[e][th], sg0_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sg0()
{
  //Get energy bin for sigma0 for given global energy
  Double_t Min = 1e38;
  Int_t e0 = 0;

  for(Int_t e=0; e<sg0_bin; e++)
    if(fabs(sg0_en[e] - gEnergy) < Min)
    {
      Min = fabs(sg0_en[e] - gEnergy);
      e0 = e;
    }
  return e0;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sg0(Double_t Energy)
{
  //Check for existing energy bin
  Int_t e0 = sg0_bin;

  for(Int_t e=0; e<sg0_bin; e++)
    if(sg0_en[e]==Energy)
      e0 = e;

  return e0;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sg0()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sg0;
  Int_t e0 = GetEnergyBin_sg0();

  //Calculate chi^2 for sigma0 data
  ChiSq_sg0 = 0.0;
  Omega = sg0_en[e0];
  for(Int_t th=0; th<sg0_pts[e0]; th++)
  {
    Theta = sg0_th[e0][th];
    Meas  = sg0_val[e0][th]*f_obs[SIG_0];
    Error = sg0_err[e0][th]*f_obs[SIG_0];
    Theo  = sigma0(Theta, Omega);
    //printf("sg0: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sg0+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sg0_wt[e0]*ChiSq_sg0; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

void Sort_sg0(Int_t l, Int_t r) //Quicksort implementation on sg0 data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sg0_en[++i] < sg0_en[r]);
     while((sg0_en[--j] > sg0_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sg0_en[i],  &sg0_en[j]);
     Swap(&sg0_wt[i],  &sg0_wt[j]);
     Swap(&sg0_sy[i],  &sg0_sy[j]);
     Swap(&sg0_pts[i], &sg0_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sg0_val[i][n], &sg0_val[j][n]);
        Swap(&sg0_err[i][n], &sg0_err[j][n]);
        Swap(&sg0_th[i][n],  &sg0_th[j][n]);
     }
   }
   Swap(&sg0_en[i],  &sg0_en[r]);
   Swap(&sg0_wt[i],  &sg0_wt[r]);
   Swap(&sg0_sy[i],  &sg0_sy[r]);
   Swap(&sg0_pts[i], &sg0_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sg0_val[i][n], &sg0_val[r][n]);
      Swap(&sg0_err[i][n], &sg0_err[r][n]);
      Swap(&sg0_th[i][n],  &sg0_th[r][n]);
   }

   Sort_sg0(l, i-1);
   Sort_sg0(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sg0()
{
  Int_t e0 = GetEnergyBin_sg0();
  Double_t Scale_sg0 = (1.0*sg0_pts[e0])*(f_obs[SIG_0]-1.0)*(f_obs[SIG_0]-1.0)/(sg0_sy[e0]*sg0_sy[e0]);
  return Scale_sg0;
}

//-----------------------------------------------------------------------------
