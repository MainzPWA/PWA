#include "Provide_sgP.h"

//-----------------------------------------------------------------------------

void Load_sgP(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaP, DsigmaP, EsigmaP;
  Char_t Ident[256];
  FILE* File_sgP;

  printf("Loading sgP  data from %s\n", Filename);
  File_sgP = fopen(Filename, "r");

  while(!feof(File_sgP))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgP, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgP, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgP(File_sgP, &Theta, &sigmaP, &DsigmaP, &EsigmaP)>=3)
    {
      sgP_val[sgP_bin][ThetaBin] = sigmaP;
      sgP_err[sgP_bin][ThetaBin] = DsigmaP;
      sgP_unc[sgP_bin][ThetaBin] = EsigmaP;
      sgP_th[sgP_bin][ThetaBin]  = Theta;
      if(DsigmaP!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgP_pts[sgP_bin] = ThetaBin;
    sgP_pre[sgP_bin] = Prelim;
    sgP_en[sgP_bin] = Energy;
    sgP_lo[sgP_bin] = EnergyLo;
    sgP_hi[sgP_bin] = EnergyHi;
    sgP_wt[sgP_bin] = Weight;
    sgP_sy[sgP_bin] = System;
    sgP_sc[sgP_bin] = Scale;
    strcpy(sgP_id[sgP_bin], Ident);

    //Increase energy bin counter
    sgP_bin++;
  }

  fclose(File_sgP);
  Sort_sgP(0, sgP_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgP_bin; t++) n+=sgP_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgP_bin; t++) if(sgP_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgP_bin);
  for(Int_t e=0; e<sgP_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgP_en[e], sgP_pts[e]);
    for(Int_t th=0; th<sgP_pts[e]; th++)
      printf("%f %f %f\n", sgP_th[e][th], sgP_val[e][th], sgP_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgP(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nP = 0;

  for(Int_t e=0; e<sgP_bin; e++)
    if((gEnergy > sgP_lo[e]) && (gEnergy < sgP_hi[e]) && (USE_PRELIMINARY || !sgP_pre[e]))
    {
      if(bins) bins[nP] = e;
      nP++;
    }

  return nP;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgP()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgP = 0.0;
  Int_t eP[EBINS];
  Int_t nP = GetEnergyBins_sgP(eP); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaP data
  for(Int_t n=0; n<nP; n++) //Process all found bins
  {
    Omega = sgP_en[eP[n]];
    for(Int_t th=0; th<sgP_pts[eP[n]]; th++) //Process all data points in current bin
    {
      Theta = sgP_th[eP[n]][th];
      Meas  = sgP_sc[eP[n]]*sgP_val[eP[n]][th]*f_obs[SIG_P];
      Error = sgP_sc[eP[n]]*sgP_err[eP[n]][th]*f_obs[SIG_P];
      Theo  = sigmaP(Theta, Omega);
      //printf("sgP: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgP+=(sgP_wt[eP[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgP;
}

//-----------------------------------------------------------------------------

void Sort_sgP(Int_t l, Int_t r) //Quicksort implementation on sgP data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgP_en[++i] < sgP_en[r]);
     while((sgP_en[--j] > sgP_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgP_lo[i],  &sgP_lo[j]);
     Swap(&sgP_en[i],  &sgP_en[j]);
     Swap(&sgP_hi[i],  &sgP_hi[j]);
     Swap(&sgP_wt[i],  &sgP_wt[j]);
     Swap(&sgP_sy[i],  &sgP_sy[j]);
     Swap(&sgP_sc[i],  &sgP_sc[j]);
     Swap(&sgP_pts[i], &sgP_pts[j]);
     Swap(&sgP_pre[i], &sgP_pre[j]);
     Swap(sgP_id[i],   sgP_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgP_val[i][n], &sgP_val[j][n]);
        Swap(&sgP_err[i][n], &sgP_err[j][n]);
        Swap(&sgP_unc[i][n], &sgP_unc[j][n]);
        Swap(&sgP_th[i][n],  &sgP_th[j][n]);
     }
   }
   Swap(&sgP_lo[i],  &sgP_lo[r]);
   Swap(&sgP_en[i],  &sgP_en[r]);
   Swap(&sgP_hi[i],  &sgP_hi[r]);
   Swap(&sgP_wt[i],  &sgP_wt[r]);
   Swap(&sgP_sy[i],  &sgP_sy[r]);
   Swap(&sgP_sc[i],  &sgP_sc[r]);
   Swap(&sgP_pts[i], &sgP_pts[r]);
   Swap(&sgP_pre[i], &sgP_pre[r]);
   Swap(sgP_id[i],   sgP_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgP_val[i][n], &sgP_val[r][n]);
      Swap(&sgP_err[i][n], &sgP_err[r][n]);
      Swap(&sgP_unc[i][n], &sgP_unc[r][n]);
      Swap(&sgP_th[i][n],  &sgP_th[r][n]);
   }

   Sort_sgP(l, i-1);
   Sort_sgP(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgP()
{
  Double_t Scale_sgP = 0.0;
  Int_t eP[EBINS];
  Int_t nP = GetEnergyBins_sgP(eP); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nP; n++) //Process all found bins
    Scale_sgP+=(f_obs[SIG_P]-1.0)*(f_obs[SIG_P]-1.0)*sgP_pts[eP[n]]/(sgP_sy[eP[n]]*sgP_sy[eP[n]]);

  return Scale_sgP;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgP()
{
  Int_t NPts_sgP = 0;
  Int_t eP[EBINS];
  Int_t nP = GetEnergyBins_sgP(eP); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nP; n++) //Process all found bins
    NPts_sgP+=sgP_pts[eP[n]];

  return NPts_sgP;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgP(FILE* File_sgP, Double_t* Theta, Double_t* sigmaP, Double_t* DsigmaP, Double_t* EsigmaP)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgP);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaP, DsigmaP, EsigmaP);
}

//-----------------------------------------------------------------------------
