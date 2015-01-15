#include "Provide_sgF.h"

//-----------------------------------------------------------------------------

void Load_sgF(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaF, DsigmaF, EsigmaF;
  Char_t Ident[256];
  FILE* File_sgF;

  printf("Loading sgF  data from %s\n", Filename);
  File_sgF = fopen(Filename, "r");

  while(!feof(File_sgF))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgF, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgF, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgF(File_sgF, &Theta, &sigmaF, &DsigmaF, &EsigmaF)>=3)
    {
      sgF_val[sgF_bin][ThetaBin] = sigmaF;
      sgF_err[sgF_bin][ThetaBin] = DsigmaF;
      sgF_unc[sgF_bin][ThetaBin] = EsigmaF;
      sgF_th[sgF_bin][ThetaBin]  = Theta;
      if(DsigmaF!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgF_pts[sgF_bin] = ThetaBin;
    sgF_pre[sgF_bin] = Prelim;
    sgF_en[sgF_bin] = Energy;
    sgF_lo[sgF_bin] = EnergyLo;
    sgF_hi[sgF_bin] = EnergyHi;
    sgF_wt[sgF_bin] = Weight;
    sgF_sy[sgF_bin] = System;
    sgF_sc[sgF_bin] = Scale;
    strcpy(sgF_id[sgF_bin], Ident);

    //Increase energy bin counter
    sgF_bin++;
  }

  fclose(File_sgF);
  Sort_sgF(0, sgF_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgF_bin; t++) n+=sgF_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgF_bin; t++) if(sgF_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgF_bin);
  for(Int_t e=0; e<sgF_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgF_en[e], sgF_pts[e]);
    for(Int_t th=0; th<sgF_pts[e]; th++)
      printf("%f %f %f\n", sgF_th[e][th], sgF_val[e][th], sgF_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgF(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nF = 0;

  for(Int_t e=0; e<sgF_bin; e++)
    if((gEnergy > sgF_lo[e]) && (gEnergy < sgF_hi[e]) && (USE_PRELIMINARY || !sgF_pre[e]))
    {
      if(bins) bins[nF] = e;
      nF++;
    }

  return nF;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgF()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgF = 0.0;
  Int_t eF[EBINS];
  Int_t nF = GetEnergyBins_sgF(eF); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaF data
  for(Int_t n=0; n<nF; n++) //Process all found bins
  {
    Omega = sgF_en[eF[n]];
    for(Int_t th=0; th<sgF_pts[eF[n]]; th++) //Process all data points in current bin
    {
      Theta = sgF_th[eF[n]][th];
      Meas  = sgF_sc[eF[n]]*sgF_val[eF[n]][th]*f_obs[SIG_F];
      Error = sgF_sc[eF[n]]*sgF_err[eF[n]][th]*f_obs[SIG_F];
      Theo  = sigmaF(Theta, Omega);
      //printf("sgF: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgF+=(sgF_wt[eF[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgF;
}

//-----------------------------------------------------------------------------

void Sort_sgF(Int_t l, Int_t r) //Quicksort implementation on sgF data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgF_en[++i] < sgF_en[r]);
     while((sgF_en[--j] > sgF_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgF_lo[i],  &sgF_lo[j]);
     Swap(&sgF_en[i],  &sgF_en[j]);
     Swap(&sgF_hi[i],  &sgF_hi[j]);
     Swap(&sgF_wt[i],  &sgF_wt[j]);
     Swap(&sgF_sy[i],  &sgF_sy[j]);
     Swap(&sgF_sc[i],  &sgF_sc[j]);
     Swap(&sgF_pts[i], &sgF_pts[j]);
     Swap(&sgF_pre[i], &sgF_pre[j]);
     Swap(sgF_id[i],   sgF_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgF_val[i][n], &sgF_val[j][n]);
        Swap(&sgF_err[i][n], &sgF_err[j][n]);
        Swap(&sgF_unc[i][n], &sgF_unc[j][n]);
        Swap(&sgF_th[i][n],  &sgF_th[j][n]);
     }
   }
   Swap(&sgF_lo[i],  &sgF_lo[r]);
   Swap(&sgF_en[i],  &sgF_en[r]);
   Swap(&sgF_hi[i],  &sgF_hi[r]);
   Swap(&sgF_wt[i],  &sgF_wt[r]);
   Swap(&sgF_sy[i],  &sgF_sy[r]);
   Swap(&sgF_sc[i],  &sgF_sc[r]);
   Swap(&sgF_pts[i], &sgF_pts[r]);
   Swap(&sgF_pre[i], &sgF_pre[r]);
   Swap(sgF_id[i],   sgF_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgF_val[i][n], &sgF_val[r][n]);
      Swap(&sgF_err[i][n], &sgF_err[r][n]);
      Swap(&sgF_unc[i][n], &sgF_unc[r][n]);
      Swap(&sgF_th[i][n],  &sgF_th[r][n]);
   }

   Sort_sgF(l, i-1);
   Sort_sgF(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgF()
{
  Double_t Scale_sgF = 0.0;
  Int_t eF[EBINS];
  Int_t nF = GetEnergyBins_sgF(eF); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nF; n++) //Process all found bins
    Scale_sgF+=(f_obs[SIG_F]-1.0)*(f_obs[SIG_F]-1.0)*sgF_pts[eF[n]]/(sgF_sy[eF[n]]*sgF_sy[eF[n]]);

  return Scale_sgF;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgF()
{
  Int_t NPts_sgF = 0;
  Int_t eF[EBINS];
  Int_t nF = GetEnergyBins_sgF(eF); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nF; n++) //Process all found bins
    NPts_sgF+=sgF_pts[eF[n]];

  return NPts_sgF;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgF(FILE* File_sgF, Double_t* Theta, Double_t* sigmaF, Double_t* DsigmaF, Double_t* EsigmaF)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgF);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaF, DsigmaF, EsigmaF);
}

//-----------------------------------------------------------------------------
