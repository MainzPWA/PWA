#include "Provide_sgLx.h"

//-----------------------------------------------------------------------------

void Load_sgLx(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaLx, DsigmaLx, EsigmaLx;
  Char_t Ident[256];
  FILE* File_sgLx;

  printf("Loading sgLx data from %s\n", Filename);
  File_sgLx = fopen(Filename, "r");

  while(!feof(File_sgLx))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgLx, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgLx, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgLx(File_sgLx, &Theta, &sigmaLx, &DsigmaLx, &EsigmaLx)>=3)
    {
      sgLx_val[sgLx_bin][ThetaBin] = sigmaLx;
      sgLx_err[sgLx_bin][ThetaBin] = DsigmaLx;
      sgLx_unc[sgLx_bin][ThetaBin] = EsigmaLx;
      sgLx_th[sgLx_bin][ThetaBin]  = Theta;
      if(DsigmaLx!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgLx_pts[sgLx_bin] = ThetaBin;
    sgLx_pre[sgLx_bin] = Prelim;
    sgLx_en[sgLx_bin] = Energy;
    sgLx_lo[sgLx_bin] = EnergyLo;
    sgLx_hi[sgLx_bin] = EnergyHi;
    sgLx_wt[sgLx_bin] = Weight;
    sgLx_sy[sgLx_bin] = System;
    sgLx_sc[sgLx_bin] = Scale;
    strcpy(sgLx_id[sgLx_bin], Ident);

    //Increase energy bin counter
    sgLx_bin++;
  }

  fclose(File_sgLx);
  Sort_sgLx(0, sgLx_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgLx_bin; t++) n+=sgLx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgLx_bin; t++) if(sgLx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgLx_bin);
  for(Int_t e=0; e<sgLx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgLx_en[e], sgLx_pts[e]);
    for(Int_t th=0; th<sgLx_pts[e]; th++)
      printf("%f %f %f\n", sgLx_th[e][th], sgLx_val[e][th], sgLx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgLx(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nLx = 0;

  for(Int_t e=0; e<sgLx_bin; e++)
    if((gEnergy > sgLx_lo[e]) && (gEnergy < sgLx_hi[e]) && (USE_PRELIMINARY || !sgLx_pre[e]))
    {
      if(bins) bins[nLx] = e;
      nLx++;
    }

  return nLx;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgLx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgLx = 0.0;
  Int_t eLx[EBINS];
  Int_t nLx = GetEnergyBins_sgLx(eLx); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaLx data
  for(Int_t n=0; n<nLx; n++) //Process all found bins
  {
    Omega = sgLx_en[eLx[n]];
    for(Int_t th=0; th<sgLx_pts[eLx[n]]; th++) //Process all data points in current bin
    {
      Theta = sgLx_th[eLx[n]][th];
      Meas  = sgLx_sc[eLx[n]]*sgLx_val[eLx[n]][th]*f_obs[SIG_LX];
      Error = sgLx_sc[eLx[n]]*sgLx_err[eLx[n]][th]*f_obs[SIG_LX];
      Theo  = sigmaLx(Theta, Omega);
      //printf("sgLx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgLx+=(sgLx_wt[eLx[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgLx;
}

//-----------------------------------------------------------------------------

void Sort_sgLx(Int_t l, Int_t r) //Quicksort implementation on sgLx data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgLx_en[++i] < sgLx_en[r]);
     while((sgLx_en[--j] > sgLx_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgLx_lo[i],  &sgLx_lo[j]);
     Swap(&sgLx_en[i],  &sgLx_en[j]);
     Swap(&sgLx_hi[i],  &sgLx_hi[j]);
     Swap(&sgLx_wt[i],  &sgLx_wt[j]);
     Swap(&sgLx_sy[i],  &sgLx_sy[j]);
     Swap(&sgLx_pts[i], &sgLx_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgLx_val[i][n], &sgLx_val[j][n]);
        Swap(&sgLx_err[i][n], &sgLx_err[j][n]);
        Swap(&sgLx_unc[i][n], &sgLx_unc[j][n]);
        Swap(&sgLx_th[i][n],  &sgLx_th[j][n]);
     }
   }
   Swap(&sgLx_lo[i],  &sgLx_lo[r]);
   Swap(&sgLx_en[i],  &sgLx_en[r]);
   Swap(&sgLx_hi[i],  &sgLx_hi[r]);
   Swap(&sgLx_wt[i],  &sgLx_wt[r]);
   Swap(&sgLx_sy[i],  &sgLx_sy[r]);
   Swap(&sgLx_pts[i], &sgLx_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgLx_val[i][n], &sgLx_val[r][n]);
      Swap(&sgLx_err[i][n], &sgLx_err[r][n]);
      Swap(&sgLx_unc[i][n], &sgLx_unc[r][n]);
      Swap(&sgLx_th[i][n],  &sgLx_th[r][n]);
   }

   Sort_sgLx(l, i-1);
   Sort_sgLx(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgLx()
{
  Double_t Scale_sgLx = 0.0;
  Int_t eLx[EBINS];
  Int_t nLx = GetEnergyBins_sgLx(eLx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nLx; n++) //Process all found bins
    Scale_sgLx+=(1.0*sgLx_pts[eLx[n]])*(f_obs[SIG_LX]-1.0)*(f_obs[SIG_LX]-1.0)/(sgLx_sy[eLx[n]]*sgLx_sy[eLx[n]]);

  return Scale_sgLx;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgLx()
{
  Int_t NPts_sgLx = 0;
  Int_t eLx[EBINS];
  Int_t nLx = GetEnergyBins_sgLx(eLx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nLx; n++) //Process all found bins
    NPts_sgLx+=sgLx_pts[eLx[n]];

  return NPts_sgLx;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgLx(FILE* File_sgLx, Double_t* Theta, Double_t* sigmaLx, Double_t* DsigmaLx, Double_t* EsigmaLx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgLx);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaLx, DsigmaLx, EsigmaLx);
}

//-----------------------------------------------------------------------------
