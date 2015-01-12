#include "Provide_sgE.h"

//-----------------------------------------------------------------------------

void Load_sgE(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaE, DsigmaE, EsigmaE;
  Char_t Ident[256];
  FILE* File_sgE;

  printf("Loading sgE  data from %s\n", Filename);
  File_sgE = fopen(Filename, "r");

  while(!feof(File_sgE))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgE, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgE, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgE(File_sgE, &Theta, &sigmaE, &DsigmaE, &EsigmaE)>=3)
    {
      sgE_val[sgE_bin][ThetaBin] = sigmaE;
      sgE_err[sgE_bin][ThetaBin] = DsigmaE;
      sgE_unc[sgE_bin][ThetaBin] = EsigmaE;
      sgE_th[sgE_bin][ThetaBin]  = Theta;
      if(DsigmaE!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgE_pts[sgE_bin] = ThetaBin;
    sgE_pre[sgE_bin] = Prelim;
    sgE_en[sgE_bin] = Energy;
    sgE_lo[sgE_bin] = EnergyLo;
    sgE_hi[sgE_bin] = EnergyHi;
    sgE_wt[sgE_bin] = Weight;
    sgE_sy[sgE_bin] = System;
    sgE_sc[sgE_bin] = Scale;
    strcpy(sgE_id[sgE_bin], Ident);

    //Increase energy bin counter
    sgE_bin++;
  }

  fclose(File_sgE);
  Sort_sgE(0, sgE_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgE_bin; t++) n+=sgE_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgE_bin; t++) if(sgE_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgE_bin);
  for(Int_t e=0; e<sgE_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgE_en[e], sgE_pts[e]);
    for(Int_t th=0; th<sgE_pts[e]; th++)
      printf("%f %f %f\n", sgE_th[e][th], sgE_val[e][th], sgE_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgE(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nE = 0;

  for(Int_t e=0; e<sgE_bin; e++)
    if((gEnergy > sgE_lo[e]) && (gEnergy < sgE_hi[e]) && (USE_PRELIMINARY || !sgE_pre[e]))
    {
      if(bins) bins[nE] = e;
      nE++;
    }

  return nE;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgE()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgE = 0.0;
  Int_t eE[EBINS];
  Int_t nE = GetEnergyBins_sgE(eE); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaE data
  for(Int_t n=0; n<nE; n++) //Process all found bins
  {
    Omega = sgE_en[eE[n]];
    for(Int_t th=0; th<sgE_pts[eE[n]]; th++) //Process all data points in current bin
    {
      Theta = sgE_th[eE[n]][th];
      Meas  = sgE_sc[eE[n]]*sgE_val[eE[n]][th]*f_obs[SIG_E];
      Error = sgE_sc[eE[n]]*sgE_err[eE[n]][th]*f_obs[SIG_E];
      Theo  = sigmaE(Theta, Omega);
      //printf("sgE: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgE+=(sgE_wt[eE[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgE;
}

//-----------------------------------------------------------------------------

void Sort_sgE(Int_t l, Int_t r) //Quicksort implementation on sgE data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgE_en[++i] < sgE_en[r]);
     while((sgE_en[--j] > sgE_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgE_lo[i],  &sgE_lo[j]);
     Swap(&sgE_en[i],  &sgE_en[j]);
     Swap(&sgE_hi[i],  &sgE_hi[j]);
     Swap(&sgE_wt[i],  &sgE_wt[j]);
     Swap(&sgE_sy[i],  &sgE_sy[j]);
     Swap(&sgE_pts[i], &sgE_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgE_val[i][n], &sgE_val[j][n]);
        Swap(&sgE_err[i][n], &sgE_err[j][n]);
        Swap(&sgE_unc[i][n], &sgE_unc[j][n]);
        Swap(&sgE_th[i][n],  &sgE_th[j][n]);
     }
   }
   Swap(&sgE_lo[i],  &sgE_lo[r]);
   Swap(&sgE_en[i],  &sgE_en[r]);
   Swap(&sgE_hi[i],  &sgE_hi[r]);
   Swap(&sgE_wt[i],  &sgE_wt[r]);
   Swap(&sgE_sy[i],  &sgE_sy[r]);
   Swap(&sgE_pts[i], &sgE_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgE_val[i][n], &sgE_val[r][n]);
      Swap(&sgE_err[i][n], &sgE_err[r][n]);
      Swap(&sgE_unc[i][n], &sgE_unc[r][n]);
      Swap(&sgE_th[i][n],  &sgE_th[r][n]);
   }

   Sort_sgE(l, i-1);
   Sort_sgE(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgE()
{
  Double_t Scale_sgE = 0.0;
  Int_t eE[EBINS];
  Int_t nE = GetEnergyBins_sgE(eE); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nE; n++) //Process all found bins
    Scale_sgE+=(f_obs[SIG_E]-1.0)*(f_obs[SIG_E]-1.0)*sgE_pts[eE[n]]/(sgE_sy[eE[n]]*sgE_sy[eE[n]]);

  return Scale_sgE;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgE()
{
  Int_t NPts_sgE = 0;
  Int_t eE[EBINS];
  Int_t nE = GetEnergyBins_sgE(eE); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nE; n++) //Process all found bins
    NPts_sgE+=sgE_pts[eE[n]];

  return NPts_sgE;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgE(FILE* File_sgE, Double_t* Theta, Double_t* sigmaE, Double_t* DsigmaE, Double_t* EsigmaE)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgE);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaE, DsigmaE, EsigmaE);
}

//-----------------------------------------------------------------------------
