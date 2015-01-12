#include "Provide_sgOx.h"

//-----------------------------------------------------------------------------

void Load_sgOx(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaOx, DsigmaOx, EsigmaOx;
  Char_t Ident[256];
  FILE* File_sgOx;

  printf("Loading sgOx data from %s\n", Filename);
  File_sgOx = fopen(Filename, "r");

  while(!feof(File_sgOx))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgOx, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgOx, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgOx(File_sgOx, &Theta, &sigmaOx, &DsigmaOx, &EsigmaOx)>=3)
    {
      sgOx_val[sgOx_bin][ThetaBin] = sigmaOx;
      sgOx_err[sgOx_bin][ThetaBin] = DsigmaOx;
      sgOx_unc[sgOx_bin][ThetaBin] = EsigmaOx;
      sgOx_th[sgOx_bin][ThetaBin]  = Theta;
      if(DsigmaOx!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgOx_pts[sgOx_bin] = ThetaBin;
    sgOx_pre[sgOx_bin] = Prelim;
    sgOx_en[sgOx_bin] = Energy;
    sgOx_lo[sgOx_bin] = EnergyLo;
    sgOx_hi[sgOx_bin] = EnergyHi;
    sgOx_wt[sgOx_bin] = Weight;
    sgOx_sy[sgOx_bin] = System;
    sgOx_sc[sgOx_bin] = Scale;
    strcpy(sgOx_id[sgOx_bin], Ident);

    //Increase energy bin counter
    sgOx_bin++;
  }

  fclose(File_sgOx);
  Sort_sgOx(0, sgOx_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgOx_bin; t++) n+=sgOx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgOx_bin; t++) if(sgOx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgOx_bin);
  for(Int_t e=0; e<sgOx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgOx_en[e], sgOx_pts[e]);
    for(Int_t th=0; th<sgOx_pts[e]; th++)
      printf("%f %f %f\n", sgOx_th[e][th], sgOx_val[e][th], sgOx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgOx(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nOx = 0;

  for(Int_t e=0; e<sgOx_bin; e++)
    if((gEnergy > sgOx_lo[e]) && (gEnergy < sgOx_hi[e]) && (USE_PRELIMINARY || !sgOx_pre[e]))
    {
      if(bins) bins[nOx] = e;
      nOx++;
    }

  return nOx;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgOx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgOx = 0.0;
  Int_t eOx[EBINS];
  Int_t nOx = GetEnergyBins_sgOx(eOx); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaOx data
  for(Int_t n=0; n<nOx; n++) //Process all found bins
  {
    Omega = sgOx_en[eOx[n]];
    for(Int_t th=0; th<sgOx_pts[eOx[n]]; th++) //Process all data points in current bin
    {
      Theta = sgOx_th[eOx[n]][th];
      Meas  = sgOx_sc[eOx[n]]*sgOx_val[eOx[n]][th]*f_obs[SIG_OX];
      Error = sgOx_sc[eOx[n]]*sgOx_err[eOx[n]][th]*f_obs[SIG_OX];
      Theo  = sigmaOx(Theta, Omega);
      //printf("sgOx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgOx+=(sgOx_wt[eOx[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgOx;
}

//-----------------------------------------------------------------------------

void Sort_sgOx(Int_t l, Int_t r) //Quicksort implementation on sgOx data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgOx_en[++i] < sgOx_en[r]);
     while((sgOx_en[--j] > sgOx_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgOx_lo[i],  &sgOx_lo[j]);
     Swap(&sgOx_en[i],  &sgOx_en[j]);
     Swap(&sgOx_hi[i],  &sgOx_hi[j]);
     Swap(&sgOx_wt[i],  &sgOx_wt[j]);
     Swap(&sgOx_sy[i],  &sgOx_sy[j]);
     Swap(&sgOx_pts[i], &sgOx_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgOx_val[i][n], &sgOx_val[j][n]);
        Swap(&sgOx_err[i][n], &sgOx_err[j][n]);
        Swap(&sgOx_unc[i][n], &sgOx_unc[j][n]);
        Swap(&sgOx_th[i][n],  &sgOx_th[j][n]);
     }
   }
   Swap(&sgOx_lo[i],  &sgOx_lo[r]);
   Swap(&sgOx_en[i],  &sgOx_en[r]);
   Swap(&sgOx_hi[i],  &sgOx_hi[r]);
   Swap(&sgOx_wt[i],  &sgOx_wt[r]);
   Swap(&sgOx_sy[i],  &sgOx_sy[r]);
   Swap(&sgOx_pts[i], &sgOx_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgOx_val[i][n], &sgOx_val[r][n]);
      Swap(&sgOx_err[i][n], &sgOx_err[r][n]);
      Swap(&sgOx_unc[i][n], &sgOx_unc[r][n]);
      Swap(&sgOx_th[i][n],  &sgOx_th[r][n]);
   }

   Sort_sgOx(l, i-1);
   Sort_sgOx(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgOx()
{
  Double_t Scale_sgOx = 0.0;
  Int_t eOx[EBINS];
  Int_t nOx = GetEnergyBins_sgOx(eOx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nOx; n++) //Process all found bins
    Scale_sgOx+=(1.0*sgOx_pts[eOx[n]])*(f_obs[SIG_OX]-1.0)*(f_obs[SIG_OX]-1.0)/(sgOx_sy[eOx[n]]*sgOx_sy[eOx[n]]);

  return Scale_sgOx;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgOx()
{
  Int_t NPts_sgOx = 0;
  Int_t eOx[EBINS];
  Int_t nOx = GetEnergyBins_sgOx(eOx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nOx; n++) //Process all found bins
    NPts_sgOx+=sgOx_pts[eOx[n]];

  return NPts_sgOx;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgOx(FILE* File_sgOx, Double_t* Theta, Double_t* sigmaOx, Double_t* DsigmaOx, Double_t* EsigmaOx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgOx);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaOx, DsigmaOx, EsigmaOx);
}

//-----------------------------------------------------------------------------
