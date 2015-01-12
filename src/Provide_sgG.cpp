#include "Provide_sgG.h"

//-----------------------------------------------------------------------------

void Load_sgG(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaG, DsigmaG, EsigmaG;
  Char_t Ident[256];
  FILE* File_sgG;

  printf("Loading sgG  data from %s\n", Filename);
  File_sgG = fopen(Filename, "r");

  while(!feof(File_sgG))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgG, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgG, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgG(File_sgG, &Theta, &sigmaG, &DsigmaG, &EsigmaG)>=3)
    {
      sgG_val[sgG_bin][ThetaBin] = sigmaG;
      sgG_err[sgG_bin][ThetaBin] = DsigmaG;
      sgG_unc[sgG_bin][ThetaBin] = EsigmaG;
      sgG_th[sgG_bin][ThetaBin]  = Theta;
      if(DsigmaG!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgG_pts[sgG_bin] = ThetaBin;
    sgG_pre[sgG_bin] = Prelim;
    sgG_en[sgG_bin] = Energy;
    sgG_lo[sgG_bin] = EnergyLo;
    sgG_hi[sgG_bin] = EnergyHi;
    sgG_wt[sgG_bin] = Weight;
    sgG_sy[sgG_bin] = System;
    sgG_sc[sgG_bin] = Scale;
    strcpy(sgG_id[sgG_bin], Ident);

    //Increase energy bin counter
    sgG_bin++;
  }

  fclose(File_sgG);
  Sort_sgG(0, sgG_bin-1);
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

Int_t GetEnergyBins_sgG(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nG = 0;

  for(Int_t e=0; e<sgG_bin; e++)
    if((gEnergy > sgG_lo[e]) && (gEnergy < sgG_hi[e]) && (USE_PRELIMINARY || !sgG_pre[e]))
    {
      if(bins) bins[nG] = e;
      nG++;
    }

  return nG;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgG()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgG = 0.0;
  Int_t eG[EBINS];
  Int_t nG = GetEnergyBins_sgG(eG); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaG data
  for(Int_t n=0; n<nG; n++) //Process all found bins
  {
    Omega = sgG_en[eG[n]];
    for(Int_t th=0; th<sgG_pts[eG[n]]; th++) //Process all data points in current bin
    {
      Theta = sgG_th[eG[n]][th];
      Meas  = sgG_sc[eG[n]]*sgG_val[eG[n]][th]*f_obs[SIG_G];
      Error = sgG_sc[eG[n]]*sgG_err[eG[n]][th]*f_obs[SIG_G];
      Theo  = sigmaG(Theta, Omega);
      //printf("sgG: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgG+=(sgG_wt[eG[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgG;
}

//-----------------------------------------------------------------------------

void Sort_sgG(Int_t l, Int_t r) //Quicksort implementation on sgG data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgG_en[++i] < sgG_en[r]);
     while((sgG_en[--j] > sgG_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgG_lo[i],  &sgG_lo[j]);
     Swap(&sgG_en[i],  &sgG_en[j]);
     Swap(&sgG_hi[i],  &sgG_hi[j]);
     Swap(&sgG_wt[i],  &sgG_wt[j]);
     Swap(&sgG_sy[i],  &sgG_sy[j]);
     Swap(&sgG_pts[i], &sgG_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgG_val[i][n], &sgG_val[j][n]);
        Swap(&sgG_err[i][n], &sgG_err[j][n]);
        Swap(&sgG_unc[i][n], &sgG_unc[j][n]);
        Swap(&sgG_th[i][n],  &sgG_th[j][n]);
     }
   }
   Swap(&sgG_lo[i],  &sgG_lo[r]);
   Swap(&sgG_en[i],  &sgG_en[r]);
   Swap(&sgG_hi[i],  &sgG_hi[r]);
   Swap(&sgG_wt[i],  &sgG_wt[r]);
   Swap(&sgG_sy[i],  &sgG_sy[r]);
   Swap(&sgG_pts[i], &sgG_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgG_val[i][n], &sgG_val[r][n]);
      Swap(&sgG_err[i][n], &sgG_err[r][n]);
      Swap(&sgG_unc[i][n], &sgG_unc[r][n]);
      Swap(&sgG_th[i][n],  &sgG_th[r][n]);
   }

   Sort_sgG(l, i-1);
   Sort_sgG(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgG()
{
  Double_t Scale_sgG = 0.0;
  Int_t eG[EBINS];
  Int_t nG = GetEnergyBins_sgG(eG); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nG; n++) //Process all found bins
    Scale_sgG+=(f_obs[SIG_G]-1.0)*(f_obs[SIG_G]-1.0)*sgG_pts[eG[n]]/(sgG_sy[eG[n]]*sgG_sy[eG[n]]);

  return Scale_sgG;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgG()
{
  Int_t NPts_sgG = 0;
  Int_t eG[EBINS];
  Int_t nG = GetEnergyBins_sgG(eG); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nG; n++) //Process all found bins
    NPts_sgG+=sgG_pts[eG[n]];

  return NPts_sgG;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgG(FILE* File_sgG, Double_t* Theta, Double_t* sigmaG, Double_t* DsigmaG, Double_t* EsigmaG)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgG);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaG, DsigmaG, EsigmaG);
}

//-----------------------------------------------------------------------------
