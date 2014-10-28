#include "Provide_sgCx.h"

//-----------------------------------------------------------------------------

void Load_sgCx(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaCx, DsigmaCx, EsigmaCx;
  Char_t Ident[256];
  FILE* File_sgCx;

  printf("Loading sgCx data from %s\n", Filename);
  File_sgCx = fopen(Filename, "r");

  while(!feof(File_sgCx))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgCx, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgCx, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgCx(File_sgCx, &Theta, &sigmaCx, &DsigmaCx, &EsigmaCx)>=3)
    {
      sgCx_val[sgCx_bin][ThetaBin] = sigmaCx;
      sgCx_err[sgCx_bin][ThetaBin] = DsigmaCx;
      sgCx_unc[sgCx_bin][ThetaBin] = EsigmaCx;
      sgCx_th[sgCx_bin][ThetaBin]  = Theta;
      if(DsigmaCx!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgCx_pts[sgCx_bin] = ThetaBin;
    sgCx_pre[sgCx_bin] = Prelim;
    sgCx_en[sgCx_bin] = Energy;
    sgCx_lo[sgCx_bin] = EnergyLo;
    sgCx_hi[sgCx_bin] = EnergyHi;
    sgCx_wt[sgCx_bin] = Weight;
    sgCx_sy[sgCx_bin] = System;
    sgCx_sc[sgCx_bin] = Scale;
    strcpy(sgCx_id[sgCx_bin], Ident);

    //Increase energy bin counter
    sgCx_bin++;
  }

  fclose(File_sgCx);
  Sort_sgCx(0, sgCx_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgCx_bin; t++) n+=sgCx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgCx_bin; t++) if(sgCx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgCx_bin);
  for(Int_t e=0; e<sgCx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgCx_en[e], sgCx_pts[e]);
    for(Int_t th=0; th<sgCx_pts[e]; th++)
      printf("%f %f %f\n", sgCx_th[e][th], sgCx_val[e][th], sgCx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgCx(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nCx = 0;

  for(Int_t e=0; e<sgCx_bin; e++)
    if((gEnergy > sgCx_lo[e]) && (gEnergy < sgCx_hi[e]) && (USE_PRELIMINARY || !sgCx_pre[e]))
    {
      if(bins) bins[nCx] = e;
      nCx++;
    }

  return nCx;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgCx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgCx = 0.0;
  Int_t eCx[EBINS];
  Int_t nCx = GetEnergyBins_sgCx(eCx); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaCx data
  for(Int_t n=0; n<nCx; n++) //Process all found bins
  {
    Omega = sgCx_en[eCx[n]];
    for(Int_t th=0; th<sgCx_pts[eCx[n]]; th++) //Process all data points in current bin
    {
      Theta = sgCx_th[eCx[n]][th];
      Meas  = sgCx_sc[eCx[n]]*sgCx_val[eCx[n]][th]*f_obs[SIG_0];
      Error = sgCx_sc[eCx[n]]*sgCx_err[eCx[n]][th]*f_obs[SIG_0];
      Theo  = sigmaCx(Theta, Omega);
      //printf("sgCx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgCx+=(sgCx_wt[eCx[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgCx;
}

//-----------------------------------------------------------------------------

void Sort_sgCx(Int_t l, Int_t r) //Quicksort implementation on sgCx data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgCx_en[++i] < sgCx_en[r]);
     while((sgCx_en[--j] > sgCx_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgCx_lo[i],  &sgCx_lo[j]);
     Swap(&sgCx_en[i],  &sgCx_en[j]);
     Swap(&sgCx_hi[i],  &sgCx_hi[j]);
     Swap(&sgCx_wt[i],  &sgCx_wt[j]);
     Swap(&sgCx_sy[i],  &sgCx_sy[j]);
     Swap(&sgCx_pts[i], &sgCx_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgCx_val[i][n], &sgCx_val[j][n]);
        Swap(&sgCx_err[i][n], &sgCx_err[j][n]);
        Swap(&sgCx_unc[i][n], &sgCx_unc[j][n]);
        Swap(&sgCx_th[i][n],  &sgCx_th[j][n]);
     }
   }
   Swap(&sgCx_lo[i],  &sgCx_lo[r]);
   Swap(&sgCx_en[i],  &sgCx_en[r]);
   Swap(&sgCx_hi[i],  &sgCx_hi[r]);
   Swap(&sgCx_wt[i],  &sgCx_wt[r]);
   Swap(&sgCx_sy[i],  &sgCx_sy[r]);
   Swap(&sgCx_pts[i], &sgCx_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgCx_val[i][n], &sgCx_val[r][n]);
      Swap(&sgCx_err[i][n], &sgCx_err[r][n]);
      Swap(&sgCx_unc[i][n], &sgCx_unc[r][n]);
      Swap(&sgCx_th[i][n],  &sgCx_th[r][n]);
   }

   Sort_sgCx(l, i-1);
   Sort_sgCx(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgCx()
{
  Double_t Scale_sgCx = 0.0;
  Int_t eCx[EBINS];
  Int_t nCx = GetEnergyBins_sgCx(eCx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nCx; n++) //Process all found bins
    Scale_sgCx+=(1.0*sgCx_pts[eCx[n]])*(f_obs[SIG_0]-1.0)*(f_obs[SIG_0]-1.0)/(sgCx_sy[eCx[n]]*sgCx_sy[eCx[n]]);

  return Scale_sgCx;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgCx()
{
  Int_t NPts_sgCx = 0;
  Int_t eCx[EBINS];
  Int_t nCx = GetEnergyBins_sgCx(eCx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nCx; n++) //Process all found bins
    NPts_sgCx+=sgCx_pts[eCx[n]];

  return NPts_sgCx;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgCx(FILE* File_sgCx, Double_t* Theta, Double_t* sigmaCx, Double_t* DsigmaCx, Double_t* EsigmaCx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgCx);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaCx, DsigmaCx, EsigmaCx);
}

//-----------------------------------------------------------------------------
