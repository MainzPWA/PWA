#include "Provide_sgTx.h"

//-----------------------------------------------------------------------------

void Load_sgTx(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaTx, DsigmaTx, EsigmaTx;
  Char_t Ident[256];
  FILE* File_sgTx;

  printf("Loading sgTx data from %s\n", Filename);
  File_sgTx = fopen(Filename, "r");

  while(!feof(File_sgTx))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgTx, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgTx, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgTx(File_sgTx, &Theta, &sigmaTx, &DsigmaTx, &EsigmaTx)>=3)
    {
      sgTx_val[sgTx_bin][ThetaBin] = sigmaTx;
      sgTx_err[sgTx_bin][ThetaBin] = DsigmaTx;
      sgTx_unc[sgTx_bin][ThetaBin] = EsigmaTx;
      sgTx_th[sgTx_bin][ThetaBin]  = Theta;
      if(DsigmaTx!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgTx_pts[sgTx_bin] = ThetaBin;
    sgTx_pre[sgTx_bin] = Prelim;
    sgTx_en[sgTx_bin] = Energy;
    sgTx_lo[sgTx_bin] = EnergyLo;
    sgTx_hi[sgTx_bin] = EnergyHi;
    sgTx_wt[sgTx_bin] = Weight;
    sgTx_sy[sgTx_bin] = System;
    sgTx_sc[sgTx_bin] = Scale;
    strcpy(sgTx_id[sgTx_bin], Ident);

    //Increase energy bin counter
    sgTx_bin++;
  }

  fclose(File_sgTx);
  Sort_sgTx(0, sgTx_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgTx_bin; t++) n+=sgTx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgTx_bin; t++) if(sgTx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgTx_bin);
  for(Int_t e=0; e<sgTx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgTx_en[e], sgTx_pts[e]);
    for(Int_t th=0; th<sgTx_pts[e]; th++)
      printf("%f %f %f\n", sgTx_th[e][th], sgTx_val[e][th], sgTx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgTx(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nTx = 0;

  for(Int_t e=0; e<sgTx_bin; e++)
    if((gEnergy > sgTx_lo[e]) && (gEnergy < sgTx_hi[e]) && (USE_PRELIMINARY || !sgTx_pre[e]))
    {
      if(bins) bins[nTx] = e;
      nTx++;
    }

  return nTx;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgTx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgTx = 0.0;
  Int_t eTx[EBINS];
  Int_t nTx = GetEnergyBins_sgTx(eTx); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaTx data
  for(Int_t n=0; n<nTx; n++) //Process all found bins
  {
    Omega = sgTx_en[eTx[n]];
    for(Int_t th=0; th<sgTx_pts[eTx[n]]; th++) //Process all data points in current bin
    {
      Theta = sgTx_th[eTx[n]][th];
      Meas  = sgTx_sc[eTx[n]]*sgTx_val[eTx[n]][th]*f_obs[SIG_TX];
      Error = sgTx_sc[eTx[n]]*sgTx_err[eTx[n]][th]*f_obs[SIG_TX];
      Theo  = sigmaTx(Theta, Omega);
      //printf("sgTx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgTx+=(sgTx_wt[eTx[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgTx;
}

//-----------------------------------------------------------------------------

void Sort_sgTx(Int_t l, Int_t r) //Quicksort implementation on sgTx data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgTx_en[++i] < sgTx_en[r]);
     while((sgTx_en[--j] > sgTx_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgTx_lo[i],  &sgTx_lo[j]);
     Swap(&sgTx_en[i],  &sgTx_en[j]);
     Swap(&sgTx_hi[i],  &sgTx_hi[j]);
     Swap(&sgTx_wt[i],  &sgTx_wt[j]);
     Swap(&sgTx_sy[i],  &sgTx_sy[j]);
     Swap(&sgTx_pts[i], &sgTx_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgTx_val[i][n], &sgTx_val[j][n]);
        Swap(&sgTx_err[i][n], &sgTx_err[j][n]);
        Swap(&sgTx_unc[i][n], &sgTx_unc[j][n]);
        Swap(&sgTx_th[i][n],  &sgTx_th[j][n]);
     }
   }
   Swap(&sgTx_lo[i],  &sgTx_lo[r]);
   Swap(&sgTx_en[i],  &sgTx_en[r]);
   Swap(&sgTx_hi[i],  &sgTx_hi[r]);
   Swap(&sgTx_wt[i],  &sgTx_wt[r]);
   Swap(&sgTx_sy[i],  &sgTx_sy[r]);
   Swap(&sgTx_pts[i], &sgTx_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgTx_val[i][n], &sgTx_val[r][n]);
      Swap(&sgTx_err[i][n], &sgTx_err[r][n]);
      Swap(&sgTx_unc[i][n], &sgTx_unc[r][n]);
      Swap(&sgTx_th[i][n],  &sgTx_th[r][n]);
   }

   Sort_sgTx(l, i-1);
   Sort_sgTx(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgTx()
{
  Double_t Scale_sgTx = 0.0;
  Int_t eTx[EBINS];
  Int_t nTx = GetEnergyBins_sgTx(eTx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nTx; n++) //Process all found bins
    Scale_sgTx+=(f_obs[SIG_TX]-1.0)*(f_obs[SIG_TX]-1.0)*sgTx_pts[eTx[n]]/(sgTx_sy[eTx[n]]*sgTx_sy[eTx[n]]);

  return Scale_sgTx;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgTx()
{
  Int_t NPts_sgTx = 0;
  Int_t eTx[EBINS];
  Int_t nTx = GetEnergyBins_sgTx(eTx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nTx; n++) //Process all found bins
    NPts_sgTx+=sgTx_pts[eTx[n]];

  return NPts_sgTx;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgTx(FILE* File_sgTx, Double_t* Theta, Double_t* sigmaTx, Double_t* DsigmaTx, Double_t* EsigmaTx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgTx);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaTx, DsigmaTx, EsigmaTx);
}

//-----------------------------------------------------------------------------
