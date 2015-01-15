#include "Provide_sgOz.h"

//-----------------------------------------------------------------------------

void Load_sgOz(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaOz, DsigmaOz, EsigmaOz;
  Char_t Ident[256];
  FILE* File_sgOz;

  printf("Loading sgOz data from %s\n", Filename);
  File_sgOz = fopen(Filename, "r");

  while(!feof(File_sgOz))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgOz, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgOz, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgOz(File_sgOz, &Theta, &sigmaOz, &DsigmaOz, &EsigmaOz)>=3)
    {
      sgOz_val[sgOz_bin][ThetaBin] = sigmaOz;
      sgOz_err[sgOz_bin][ThetaBin] = DsigmaOz;
      sgOz_unc[sgOz_bin][ThetaBin] = EsigmaOz;
      sgOz_th[sgOz_bin][ThetaBin]  = Theta;
      if(DsigmaOz!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgOz_pts[sgOz_bin] = ThetaBin;
    sgOz_pre[sgOz_bin] = Prelim;
    sgOz_en[sgOz_bin] = Energy;
    sgOz_lo[sgOz_bin] = EnergyLo;
    sgOz_hi[sgOz_bin] = EnergyHi;
    sgOz_wt[sgOz_bin] = Weight;
    sgOz_sy[sgOz_bin] = System;
    sgOz_sc[sgOz_bin] = Scale;
    strcpy(sgOz_id[sgOz_bin], Ident);

    //Increase energy bin counter
    sgOz_bin++;
  }

  fclose(File_sgOz);
  Sort_sgOz(0, sgOz_bin-1);
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

Int_t GetEnergyBins_sgOz(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nOz = 0;

  for(Int_t e=0; e<sgOz_bin; e++)
    if((gEnergy > sgOz_lo[e]) && (gEnergy < sgOz_hi[e]) && (USE_PRELIMINARY || !sgOz_pre[e]))
    {
      if(bins) bins[nOz] = e;
      nOz++;
    }

  return nOz;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgOz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgOz = 0.0;
  Int_t eOz[EBINS];
  Int_t nOz = GetEnergyBins_sgOz(eOz); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaOz data
  for(Int_t n=0; n<nOz; n++) //Process all found bins
  {
    Omega = sgOz_en[eOz[n]];
    for(Int_t th=0; th<sgOz_pts[eOz[n]]; th++) //Process all data points in current bin
    {
      Theta = sgOz_th[eOz[n]][th];
      Meas  = sgOz_sc[eOz[n]]*sgOz_val[eOz[n]][th]*f_obs[SIG_OZ];
      Error = sgOz_sc[eOz[n]]*sgOz_err[eOz[n]][th]*f_obs[SIG_OZ];
      Theo  = sigmaOz(Theta, Omega);
      //printf("sgOz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgOz+=(sgOz_wt[eOz[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgOz;
}

//-----------------------------------------------------------------------------

void Sort_sgOz(Int_t l, Int_t r) //Quicksort implementation on sgOz data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgOz_en[++i] < sgOz_en[r]);
     while((sgOz_en[--j] > sgOz_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgOz_lo[i],  &sgOz_lo[j]);
     Swap(&sgOz_en[i],  &sgOz_en[j]);
     Swap(&sgOz_hi[i],  &sgOz_hi[j]);
     Swap(&sgOz_wt[i],  &sgOz_wt[j]);
     Swap(&sgOz_sy[i],  &sgOz_sy[j]);
     Swap(&sgOz_sc[i],  &sgOz_sc[j]);
     Swap(&sgOz_pts[i], &sgOz_pts[j]);
     Swap(&sgOz_pre[i], &sgOz_pre[j]);
     Swap(sgOz_id[i],   sgOz_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgOz_val[i][n], &sgOz_val[j][n]);
        Swap(&sgOz_err[i][n], &sgOz_err[j][n]);
        Swap(&sgOz_unc[i][n], &sgOz_unc[j][n]);
        Swap(&sgOz_th[i][n],  &sgOz_th[j][n]);
     }
   }
   Swap(&sgOz_lo[i],  &sgOz_lo[r]);
   Swap(&sgOz_en[i],  &sgOz_en[r]);
   Swap(&sgOz_hi[i],  &sgOz_hi[r]);
   Swap(&sgOz_wt[i],  &sgOz_wt[r]);
   Swap(&sgOz_sy[i],  &sgOz_sy[r]);
   Swap(&sgOz_sc[i],  &sgOz_sc[r]);
   Swap(&sgOz_pts[i], &sgOz_pts[r]);
   Swap(&sgOz_pre[i], &sgOz_pre[r]);
   Swap(sgOz_id[i],   sgOz_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgOz_val[i][n], &sgOz_val[r][n]);
      Swap(&sgOz_err[i][n], &sgOz_err[r][n]);
      Swap(&sgOz_unc[i][n], &sgOz_unc[r][n]);
      Swap(&sgOz_th[i][n],  &sgOz_th[r][n]);
   }

   Sort_sgOz(l, i-1);
   Sort_sgOz(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgOz()
{
  Double_t Scale_sgOz = 0.0;
  Int_t eOz[EBINS];
  Int_t nOz = GetEnergyBins_sgOz(eOz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nOz; n++) //Process all found bins
    Scale_sgOz+=(f_obs[SIG_OZ]-1.0)*(f_obs[SIG_OZ]-1.0)*sgOz_pts[eOz[n]]/(sgOz_sy[eOz[n]]*sgOz_sy[eOz[n]]);

  return Scale_sgOz;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgOz()
{
  Int_t NPts_sgOz = 0;
  Int_t eOz[EBINS];
  Int_t nOz = GetEnergyBins_sgOz(eOz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nOz; n++) //Process all found bins
    NPts_sgOz+=sgOz_pts[eOz[n]];

  return NPts_sgOz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgOz(FILE* File_sgOz, Double_t* Theta, Double_t* sigmaOz, Double_t* DsigmaOz, Double_t* EsigmaOz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgOz);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaOz, DsigmaOz, EsigmaOz);
}

//-----------------------------------------------------------------------------
