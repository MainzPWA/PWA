#include "Provide_sgLz.h"

//-----------------------------------------------------------------------------

void Load_sgLz(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaLz, DsigmaLz, EsigmaLz;
  Char_t Ident[256];
  FILE* File_sgLz;

  printf("Loading sgLz data from %s\n", Filename);
  File_sgLz = fopen(Filename, "r");

  while(!feof(File_sgLz))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgLz, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgLz, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgLz(File_sgLz, &Theta, &sigmaLz, &DsigmaLz, &EsigmaLz)>=3)
    {
      sgLz_val[sgLz_bin][ThetaBin] = sigmaLz;
      sgLz_err[sgLz_bin][ThetaBin] = DsigmaLz;
      sgLz_unc[sgLz_bin][ThetaBin] = EsigmaLz;
      sgLz_th[sgLz_bin][ThetaBin]  = Theta;
      if(DsigmaLz!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgLz_pts[sgLz_bin] = ThetaBin;
    sgLz_pre[sgLz_bin] = Prelim;
    sgLz_en[sgLz_bin] = Energy;
    sgLz_lo[sgLz_bin] = EnergyLo;
    sgLz_hi[sgLz_bin] = EnergyHi;
    sgLz_wt[sgLz_bin] = Weight;
    sgLz_sy[sgLz_bin] = System;
    sgLz_sc[sgLz_bin] = Scale;
    strcpy(sgLz_id[sgLz_bin], Ident);

    //Increase energy bin counter
    sgLz_bin++;
  }

  fclose(File_sgLz);
  Sort_sgLz(0, sgLz_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgLz_bin; t++) n+=sgLz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgLz_bin; t++) if(sgLz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgLz_bin);
  for(Int_t e=0; e<sgLz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgLz_en[e], sgLz_pts[e]);
    for(Int_t th=0; th<sgLz_pts[e]; th++)
      printf("%f %f %f\n", sgLz_th[e][th], sgLz_val[e][th], sgLz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgLz(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nLz = 0;

  for(Int_t e=0; e<sgLz_bin; e++)
    if((gEnergy > sgLz_lo[e]) && (gEnergy < sgLz_hi[e]) && (USE_PRELIMINARY || !sgLz_pre[e]))
    {
      if(bins) bins[nLz] = e;
      nLz++;
    }

  return nLz;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgLz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgLz = 0.0;
  Int_t eLz[EBINS];
  Int_t nLz = GetEnergyBins_sgLz(eLz); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaLz data
  for(Int_t n=0; n<nLz; n++) //Process all found bins
  {
    Omega = sgLz_en[eLz[n]];
    for(Int_t th=0; th<sgLz_pts[eLz[n]]; th++) //Process all data points in current bin
    {
      Theta = sgLz_th[eLz[n]][th];
      Meas  = sgLz_sc[eLz[n]]*sgLz_val[eLz[n]][th]*f_obs[SIG_LZ];
      Error = sgLz_sc[eLz[n]]*sgLz_err[eLz[n]][th]*f_obs[SIG_LZ];
      Theo  = sigmaLz(Theta, Omega);
      //printf("sgLz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgLz+=(sgLz_wt[eLz[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgLz;
}

//-----------------------------------------------------------------------------

void Sort_sgLz(Int_t l, Int_t r) //Quicksort implementation on sgLz data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgLz_en[++i] < sgLz_en[r]);
     while((sgLz_en[--j] > sgLz_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgLz_lo[i],  &sgLz_lo[j]);
     Swap(&sgLz_en[i],  &sgLz_en[j]);
     Swap(&sgLz_hi[i],  &sgLz_hi[j]);
     Swap(&sgLz_wt[i],  &sgLz_wt[j]);
     Swap(&sgLz_sy[i],  &sgLz_sy[j]);
     Swap(&sgLz_pts[i], &sgLz_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgLz_val[i][n], &sgLz_val[j][n]);
        Swap(&sgLz_err[i][n], &sgLz_err[j][n]);
        Swap(&sgLz_unc[i][n], &sgLz_unc[j][n]);
        Swap(&sgLz_th[i][n],  &sgLz_th[j][n]);
     }
   }
   Swap(&sgLz_lo[i],  &sgLz_lo[r]);
   Swap(&sgLz_en[i],  &sgLz_en[r]);
   Swap(&sgLz_hi[i],  &sgLz_hi[r]);
   Swap(&sgLz_wt[i],  &sgLz_wt[r]);
   Swap(&sgLz_sy[i],  &sgLz_sy[r]);
   Swap(&sgLz_pts[i], &sgLz_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgLz_val[i][n], &sgLz_val[r][n]);
      Swap(&sgLz_err[i][n], &sgLz_err[r][n]);
      Swap(&sgLz_unc[i][n], &sgLz_unc[r][n]);
      Swap(&sgLz_th[i][n],  &sgLz_th[r][n]);
   }

   Sort_sgLz(l, i-1);
   Sort_sgLz(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgLz()
{
  Double_t Scale_sgLz = 0.0;
  Int_t eLz[EBINS];
  Int_t nLz = GetEnergyBins_sgLz(eLz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nLz; n++) //Process all found bins
    Scale_sgLz+=(f_obs[SIG_LZ]-1.0)*(f_obs[SIG_LZ]-1.0)*sgLz_pts[eLz[n]]/(sgLz_sy[eLz[n]]*sgLz_sy[eLz[n]]);

  return Scale_sgLz;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgLz()
{
  Int_t NPts_sgLz = 0;
  Int_t eLz[EBINS];
  Int_t nLz = GetEnergyBins_sgLz(eLz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nLz; n++) //Process all found bins
    NPts_sgLz+=sgLz_pts[eLz[n]];

  return NPts_sgLz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgLz(FILE* File_sgLz, Double_t* Theta, Double_t* sigmaLz, Double_t* DsigmaLz, Double_t* EsigmaLz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgLz);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaLz, DsigmaLz, EsigmaLz);
}

//-----------------------------------------------------------------------------
