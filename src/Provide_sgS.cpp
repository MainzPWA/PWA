#include "Provide_sgS.h"

//-----------------------------------------------------------------------------

void Load_sgS(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaS, DsigmaS, EsigmaS;
  Char_t Ident[256];
  FILE* File_sgS;

  printf("Loading sgS  data from %s\n", Filename);
  File_sgS = fopen(Filename, "r");

  while(!feof(File_sgS))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgS, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgS, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgS(File_sgS, &Theta, &sigmaS, &DsigmaS, &EsigmaS)>=3)
    {
      sgS_val[sgS_bin][ThetaBin] = sigmaS;
      sgS_err[sgS_bin][ThetaBin] = DsigmaS;
      sgS_unc[sgS_bin][ThetaBin] = EsigmaS;
      sgS_th[sgS_bin][ThetaBin]  = Theta;
      if(DsigmaS!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgS_pts[sgS_bin] = ThetaBin;
    sgS_pre[sgS_bin] = Prelim;
    sgS_en[sgS_bin] = Energy;
    sgS_lo[sgS_bin] = EnergyLo;
    sgS_hi[sgS_bin] = EnergyHi;
    sgS_wt[sgS_bin] = Weight;
    sgS_sy[sgS_bin] = System;
    sgS_sc[sgS_bin] = Scale;
    strcpy(sgS_id[sgS_bin], Ident);

    //Increase energy bin counter
    sgS_bin++;
  }

  fclose(File_sgS);
  Sort_sgS(0, sgS_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgS_bin; t++) n+=sgS_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgS_bin; t++) if(sgS_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgS_bin);
  for(Int_t e=0; e<sgS_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgS_en[e], sgS_pts[e]);
    for(Int_t th=0; th<sgS_pts[e]; th++)
      printf("%f %f %f\n", sgS_th[e][th], sgS_val[e][th], sgS_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgS(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nS = 0;

  for(Int_t e=0; e<sgS_bin; e++)
    if((gEnergy > sgS_lo[e]) && (gEnergy < sgS_hi[e]) && (USE_PRELIMINARY || !sgS_pre[e]))
    {
      if(bins) bins[nS] = e;
      nS++;
    }

  return nS;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgS()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgS = 0.0;
  Int_t eS[EBINS];
  Int_t nS = GetEnergyBins_sgS(eS); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaS data
  for(Int_t n=0; n<nS; n++) //Process all found bins
  {
    Omega = sgS_en[eS[n]];
    for(Int_t th=0; th<sgS_pts[eS[n]]; th++) //Process all data points in current bin
    {
      Theta = sgS_th[eS[n]][th];
      Meas  = sgS_sc[eS[n]]*sgS_val[eS[n]][th]*f_obs[SIG_S];
      Error = sgS_sc[eS[n]]*sgS_err[eS[n]][th]*f_obs[SIG_S];
      Theo  = sigmaS(Theta, Omega);
      //printf("sgS: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgS+=(sgS_wt[eS[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgS;
}

//-----------------------------------------------------------------------------

void Sort_sgS(Int_t l, Int_t r) //Quicksort implementation on sgS data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgS_en[++i] < sgS_en[r]);
     while((sgS_en[--j] > sgS_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgS_lo[i],  &sgS_lo[j]);
     Swap(&sgS_en[i],  &sgS_en[j]);
     Swap(&sgS_hi[i],  &sgS_hi[j]);
     Swap(&sgS_wt[i],  &sgS_wt[j]);
     Swap(&sgS_sy[i],  &sgS_sy[j]);
     Swap(&sgS_pts[i], &sgS_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgS_val[i][n], &sgS_val[j][n]);
        Swap(&sgS_err[i][n], &sgS_err[j][n]);
        Swap(&sgS_unc[i][n], &sgS_unc[j][n]);
        Swap(&sgS_th[i][n],  &sgS_th[j][n]);
     }
   }
   Swap(&sgS_lo[i],  &sgS_lo[r]);
   Swap(&sgS_en[i],  &sgS_en[r]);
   Swap(&sgS_hi[i],  &sgS_hi[r]);
   Swap(&sgS_wt[i],  &sgS_wt[r]);
   Swap(&sgS_sy[i],  &sgS_sy[r]);
   Swap(&sgS_pts[i], &sgS_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgS_val[i][n], &sgS_val[r][n]);
      Swap(&sgS_err[i][n], &sgS_err[r][n]);
      Swap(&sgS_unc[i][n], &sgS_unc[r][n]);
      Swap(&sgS_th[i][n],  &sgS_th[r][n]);
   }

   Sort_sgS(l, i-1);
   Sort_sgS(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgS()
{
  Double_t Scale_sgS = 0.0;
  Int_t eS[EBINS];
  Int_t nS = GetEnergyBins_sgS(eS); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nS; n++) //Process all found bins
    Scale_sgS+=(1.0*sgS_pts[eS[n]])*(f_obs[SIG_S]-1.0)*(f_obs[SIG_S]-1.0)/(sgS_sy[eS[n]]*sgS_sy[eS[n]]);

  return Scale_sgS;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgS()
{
  Int_t NPts_sgS = 0;
  Int_t eS[EBINS];
  Int_t nS = GetEnergyBins_sgS(eS); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nS; n++) //Process all found bins
    NPts_sgS+=sgS_pts[eS[n]];

  return NPts_sgS;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgS(FILE* File_sgS, Double_t* Theta, Double_t* sigmaS, Double_t* DsigmaS, Double_t* EsigmaS)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgS);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaS, DsigmaS, EsigmaS);
}

//-----------------------------------------------------------------------------
