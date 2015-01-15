#include "Provide_sgTz.h"

//-----------------------------------------------------------------------------

void Load_sgTz(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaTz, DsigmaTz, EsigmaTz;
  Char_t Ident[256];
  FILE* File_sgTz;

  printf("Loading sgTz data from %s\n", Filename);
  File_sgTz = fopen(Filename, "r");

  while(!feof(File_sgTz))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgTz, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgTz, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgTz(File_sgTz, &Theta, &sigmaTz, &DsigmaTz, &EsigmaTz)>=3)
    {
      sgTz_val[sgTz_bin][ThetaBin] = sigmaTz;
      sgTz_err[sgTz_bin][ThetaBin] = DsigmaTz;
      sgTz_unc[sgTz_bin][ThetaBin] = EsigmaTz;
      sgTz_th[sgTz_bin][ThetaBin]  = Theta;
      if(DsigmaTz!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgTz_pts[sgTz_bin] = ThetaBin;
    sgTz_pre[sgTz_bin] = Prelim;
    sgTz_en[sgTz_bin] = Energy;
    sgTz_lo[sgTz_bin] = EnergyLo;
    sgTz_hi[sgTz_bin] = EnergyHi;
    sgTz_wt[sgTz_bin] = Weight;
    sgTz_sy[sgTz_bin] = System;
    sgTz_sc[sgTz_bin] = Scale;
    strcpy(sgTz_id[sgTz_bin], Ident);

    //Increase energy bin counter
    sgTz_bin++;
  }

  fclose(File_sgTz);
  Sort_sgTz(0, sgTz_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgTz_bin; t++) n+=sgTz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgTz_bin; t++) if(sgTz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgTz_bin);
  for(Int_t e=0; e<sgTz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgTz_en[e], sgTz_pts[e]);
    for(Int_t th=0; th<sgTz_pts[e]; th++)
      printf("%f %f %f\n", sgTz_th[e][th], sgTz_val[e][th], sgTz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgTz(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nTz = 0;

  for(Int_t e=0; e<sgTz_bin; e++)
    if((gEnergy > sgTz_lo[e]) && (gEnergy < sgTz_hi[e]) && (USE_PRELIMINARY || !sgTz_pre[e]))
    {
      if(bins) bins[nTz] = e;
      nTz++;
    }

  return nTz;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgTz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgTz = 0.0;
  Int_t eTz[EBINS];
  Int_t nTz = GetEnergyBins_sgTz(eTz); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaTz data
  for(Int_t n=0; n<nTz; n++) //Process all found bins
  {
    Omega = sgTz_en[eTz[n]];
    for(Int_t th=0; th<sgTz_pts[eTz[n]]; th++) //Process all data points in current bin
    {
      Theta = sgTz_th[eTz[n]][th];
      Meas  = sgTz_sc[eTz[n]]*sgTz_val[eTz[n]][th]*f_obs[SIG_TZ];
      Error = sgTz_sc[eTz[n]]*sgTz_err[eTz[n]][th]*f_obs[SIG_TZ];
      Theo  = sigmaTz(Theta, Omega);
      //printf("sgTz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgTz+=(sgTz_wt[eTz[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgTz;
}

//-----------------------------------------------------------------------------

void Sort_sgTz(Int_t l, Int_t r) //Quicksort implementation on sgTz data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgTz_en[++i] < sgTz_en[r]);
     while((sgTz_en[--j] > sgTz_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgTz_lo[i],  &sgTz_lo[j]);
     Swap(&sgTz_en[i],  &sgTz_en[j]);
     Swap(&sgTz_hi[i],  &sgTz_hi[j]);
     Swap(&sgTz_wt[i],  &sgTz_wt[j]);
     Swap(&sgTz_sy[i],  &sgTz_sy[j]);
     Swap(&sgTz_sc[i],  &sgTz_sc[j]);
     Swap(&sgTz_pts[i], &sgTz_pts[j]);
     Swap(&sgTz_pre[i], &sgTz_pre[j]);
     Swap(sgTz_id[i],   sgTz_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgTz_val[i][n], &sgTz_val[j][n]);
        Swap(&sgTz_err[i][n], &sgTz_err[j][n]);
        Swap(&sgTz_unc[i][n], &sgTz_unc[j][n]);
        Swap(&sgTz_th[i][n],  &sgTz_th[j][n]);
     }
   }
   Swap(&sgTz_lo[i],  &sgTz_lo[r]);
   Swap(&sgTz_en[i],  &sgTz_en[r]);
   Swap(&sgTz_hi[i],  &sgTz_hi[r]);
   Swap(&sgTz_wt[i],  &sgTz_wt[r]);
   Swap(&sgTz_sy[i],  &sgTz_sy[r]);
   Swap(&sgTz_sc[i],  &sgTz_sc[r]);
   Swap(&sgTz_pts[i], &sgTz_pts[r]);
   Swap(&sgTz_pre[i], &sgTz_pre[r]);
   Swap(sgTz_id[i],   sgTz_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgTz_val[i][n], &sgTz_val[r][n]);
      Swap(&sgTz_err[i][n], &sgTz_err[r][n]);
      Swap(&sgTz_unc[i][n], &sgTz_unc[r][n]);
      Swap(&sgTz_th[i][n],  &sgTz_th[r][n]);
   }

   Sort_sgTz(l, i-1);
   Sort_sgTz(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgTz()
{
  Double_t Scale_sgTz = 0.0;
  Int_t eTz[EBINS];
  Int_t nTz = GetEnergyBins_sgTz(eTz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nTz; n++) //Process all found bins
    Scale_sgTz+=(f_obs[SIG_TZ]-1.0)*(f_obs[SIG_TZ]-1.0)*sgTz_pts[eTz[n]]/(sgTz_sy[eTz[n]]*sgTz_sy[eTz[n]]);

  return Scale_sgTz;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgTz()
{
  Int_t NPts_sgTz = 0;
  Int_t eTz[EBINS];
  Int_t nTz = GetEnergyBins_sgTz(eTz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nTz; n++) //Process all found bins
    NPts_sgTz+=sgTz_pts[eTz[n]];

  return NPts_sgTz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgTz(FILE* File_sgTz, Double_t* Theta, Double_t* sigmaTz, Double_t* DsigmaTz, Double_t* EsigmaTz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgTz);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaTz, DsigmaTz, EsigmaTz);
}

//-----------------------------------------------------------------------------
