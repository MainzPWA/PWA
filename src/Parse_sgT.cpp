#include "Parse_sgT.h"

//-----------------------------------------------------------------------------

void Parse_sgT()
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, sigmaT, DsigmaT, Dummy;
  FILE* File_sgT;

  printf("Loading sgT  data... ");
  File_sgT = fopen("data/sgT.txt", "r");

  sgT_bin = 0;
  while(!feof(File_sgT))
  {
    //Get beam energy
    if(fscanf(File_sgT, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_sgT(Energy);
    if(EnergyBin==sgT_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = sgT_pts[EnergyBin]; //..hence we append to the existing energy bin

    while(fscanf(File_sgT, "%lf %lf %lf\n", &Theta, &sigmaT, &DsigmaT, &Dummy)>=3)
    {
      sgT_val[EnergyBin][ThetaBin] = sigmaT;
      sgT_err[EnergyBin][ThetaBin] = DsigmaT;
      sgT_th[EnergyBin][ThetaBin]  = Theta;
      if(DsigmaT!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgT);

    sgT_pts[EnergyBin] = ThetaBin;
    sgT_en[EnergyBin] = Energy;
    sgT_wt[EnergyBin] = Weight;
    sgT_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==sgT_bin)
      sgT_bin++;
  }

  fclose(File_sgT);
  Sort_sgT(0, sgT_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgT_bin; t++) n+=sgT_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgT_bin; t++) if(sgT_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgT_bin);
  for(Int_t e=0; e<sgT_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgT_en[e], sgT_pts[e]);
    for(Int_t th=0; th<sgT_pts[e]; th++)
      printf("%f %f %f\n", sgT_th[e][th], sgT_val[e][th], sgT_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sgT()
{
  //Get energy bin for sigma0*T for given global energy
  Double_t Min = 1e38;
  Int_t eT = 0;

  for(Int_t e=0; e<sgT_bin; e++)
    if(fabs(sgT_en[e] - gEnergy) < Min)
    {
      Min = fabs(sgT_en[e] - gEnergy);
      eT = e;
    }
  return eT;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_sgT(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eT = sgT_bin;

  for(Int_t e=0; e<sgT_bin; e++)
    if(sgT_en[e]==Energy)
      eT = e;

  return eT;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgT()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgT;
  Int_t eT = GetEnergyBin_sgT();

  //Calculate chi^2 for sigma0*T data
  ChiSq_sgT = 0.0;
  Omega = sgT_en[eT];
  for(Int_t th=0; th<sgT_pts[eT]; th++)
  {
    Theta = sgT_th[eT][th];
    Meas  = sgT_val[eT][th]*f_obs[SIG_T];
    Error = sgT_err[eT][th]*f_obs[SIG_T];
    Theo  = sigmaT(Theta, Omega);
    //printf("sgT: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_sgT+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return sgT_wt[eT]*ChiSq_sgT; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

void Sort_sgT(Int_t l, Int_t r) //Quicksort implementation on sgT data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgT_en[++i] < sgT_en[r]);
     while((sgT_en[--j] > sgT_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgT_en[i],  &sgT_en[j]);
     Swap(&sgT_wt[i],  &sgT_wt[j]);
     Swap(&sgT_sy[i],  &sgT_sy[j]);
     Swap(&sgT_pts[i], &sgT_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgT_val[i][n], &sgT_val[j][n]);
        Swap(&sgT_err[i][n], &sgT_err[j][n]);
        Swap(&sgT_th[i][n],  &sgT_th[j][n]);
     }
   }
   Swap(&sgT_en[i],  &sgT_en[r]);
   Swap(&sgT_wt[i],  &sgT_wt[r]);
   Swap(&sgT_sy[i],  &sgT_sy[r]);
   Swap(&sgT_pts[i], &sgT_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgT_val[i][n], &sgT_val[r][n]);
      Swap(&sgT_err[i][n], &sgT_err[r][n]);
      Swap(&sgT_th[i][n],  &sgT_th[r][n]);
   }

   Sort_sgT(l, i-1);
   Sort_sgT(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgT()
{
  Int_t eT = GetEnergyBin_sgT();
  Double_t Scale_sgT = (1.0*sgT_pts[eT])*(f_obs[SIG_T]-1.0)*(f_obs[SIG_T]-1.0)/(sgT_sy[eT]*sgT_sy[eT]);
  return Scale_sgT;
}

//-----------------------------------------------------------------------------
