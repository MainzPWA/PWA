#include "Provide_sgCz.h"

//-----------------------------------------------------------------------------

void Load_sgCz(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaCz, DsigmaCz, EsigmaCz;
  Char_t Ident[256];
  FILE* File_sgCz;

  printf("Loading sgCz data from %s\n", Filename);
  File_sgCz = fopen(Filename, "r");

  while(!feof(File_sgCz))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgCz, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgCz, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgCz(File_sgCz, &Theta, &sigmaCz, &DsigmaCz, &EsigmaCz)>=3)
    {
      sgCz_val[sgCz_bin][ThetaBin] = sigmaCz;
      sgCz_err[sgCz_bin][ThetaBin] = DsigmaCz;
      sgCz_unc[sgCz_bin][ThetaBin] = EsigmaCz;
      sgCz_th[sgCz_bin][ThetaBin]  = Theta;
      if(DsigmaCz!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgCz_pts[sgCz_bin] = ThetaBin;
    sgCz_pre[sgCz_bin] = Prelim;
    sgCz_en[sgCz_bin] = Energy;
    sgCz_lo[sgCz_bin] = EnergyLo;
    sgCz_hi[sgCz_bin] = EnergyHi;
    sgCz_wt[sgCz_bin] = Weight;
    sgCz_sy[sgCz_bin] = System;
    sgCz_sc[sgCz_bin] = Scale;
    strcpy(sgCz_id[sgCz_bin], Ident);

    //Increase energy bin counter
    sgCz_bin++;
  }

  fclose(File_sgCz);
  Sort_sgCz(0, sgCz_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgCz_bin; t++) n+=sgCz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgCz_bin; t++) if(sgCz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgCz_bin);
  for(Int_t e=0; e<sgCz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgCz_en[e], sgCz_pts[e]);
    for(Int_t th=0; th<sgCz_pts[e]; th++)
      printf("%f %f %f\n", sgCz_th[e][th], sgCz_val[e][th], sgCz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgCz(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nCz = 0;

  for(Int_t e=0; e<sgCz_bin; e++)
    if((gEnergy > sgCz_lo[e]) && (gEnergy < sgCz_hi[e]) && (USE_PRELIMINARY || !sgCz_pre[e]))
    {
      if(bins) bins[nCz] = e;
      nCz++;
    }

  return nCz;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgCz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgCz = 0.0;
  Int_t eCz[EBINS];
  Int_t nCz = GetEnergyBins_sgCz(eCz); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaCz data
  for(Int_t n=0; n<nCz; n++) //Process all found bins
  {
    Omega = sgCz_en[eCz[n]];
    for(Int_t th=0; th<sgCz_pts[eCz[n]]; th++) //Process all data points in current bin
    {
      Theta = sgCz_th[eCz[n]][th];
      Meas  = sgCz_sc[eCz[n]]*sgCz_val[eCz[n]][th]*f_obs[SIG_0];
      Error = sgCz_sc[eCz[n]]*sgCz_err[eCz[n]][th]*f_obs[SIG_0];
      Theo  = sigmaCz(Theta, Omega);
      //printf("sgCz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgCz+=(sgCz_wt[eCz[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgCz;
}

//-----------------------------------------------------------------------------

void Sort_sgCz(Int_t l, Int_t r) //Quicksort implementation on sgCz data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgCz_en[++i] < sgCz_en[r]);
     while((sgCz_en[--j] > sgCz_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgCz_lo[i],  &sgCz_lo[j]);
     Swap(&sgCz_en[i],  &sgCz_en[j]);
     Swap(&sgCz_hi[i],  &sgCz_hi[j]);
     Swap(&sgCz_wt[i],  &sgCz_wt[j]);
     Swap(&sgCz_sy[i],  &sgCz_sy[j]);
     Swap(&sgCz_pts[i], &sgCz_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgCz_val[i][n], &sgCz_val[j][n]);
        Swap(&sgCz_err[i][n], &sgCz_err[j][n]);
        Swap(&sgCz_unc[i][n], &sgCz_unc[j][n]);
        Swap(&sgCz_th[i][n],  &sgCz_th[j][n]);
     }
   }
   Swap(&sgCz_lo[i],  &sgCz_lo[r]);
   Swap(&sgCz_en[i],  &sgCz_en[r]);
   Swap(&sgCz_hi[i],  &sgCz_hi[r]);
   Swap(&sgCz_wt[i],  &sgCz_wt[r]);
   Swap(&sgCz_sy[i],  &sgCz_sy[r]);
   Swap(&sgCz_pts[i], &sgCz_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgCz_val[i][n], &sgCz_val[r][n]);
      Swap(&sgCz_err[i][n], &sgCz_err[r][n]);
      Swap(&sgCz_unc[i][n], &sgCz_unc[r][n]);
      Swap(&sgCz_th[i][n],  &sgCz_th[r][n]);
   }

   Sort_sgCz(l, i-1);
   Sort_sgCz(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgCz()
{
  Double_t Scale_sgCz = 0.0;
  Int_t eCz[EBINS];
  Int_t nCz = GetEnergyBins_sgCz(eCz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nCz; n++) //Process all found bins
    Scale_sgCz+=(1.0*sgCz_pts[eCz[n]])*(f_obs[SIG_0]-1.0)*(f_obs[SIG_0]-1.0)/(sgCz_sy[eCz[n]]*sgCz_sy[eCz[n]]);

  return Scale_sgCz;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgCz()
{
  Int_t NPts_sgCz = 0;
  Int_t eCz[EBINS];
  Int_t nCz = GetEnergyBins_sgCz(eCz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nCz; n++) //Process all found bins
    NPts_sgCz+=sgCz_pts[eCz[n]];

  return NPts_sgCz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgCz(FILE* File_sgCz, Double_t* Theta, Double_t* sigmaCz, Double_t* DsigmaCz, Double_t* EsigmaCz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgCz);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaCz, DsigmaCz, EsigmaCz);
}

//-----------------------------------------------------------------------------
