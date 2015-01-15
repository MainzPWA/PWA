#include "Provide_sgH.h"

//-----------------------------------------------------------------------------

void Load_sgH(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaH, DsigmaH, EsigmaH;
  Char_t Ident[256];
  FILE* File_sgH;

  printf("Loading sgH  data from %s\n", Filename);
  File_sgH = fopen(Filename, "r");

  while(!feof(File_sgH))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgH, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgH, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgH(File_sgH, &Theta, &sigmaH, &DsigmaH, &EsigmaH)>=3)
    {
      sgH_val[sgH_bin][ThetaBin] = sigmaH;
      sgH_err[sgH_bin][ThetaBin] = DsigmaH;
      sgH_unc[sgH_bin][ThetaBin] = EsigmaH;
      sgH_th[sgH_bin][ThetaBin]  = Theta;
      if(DsigmaH!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgH_pts[sgH_bin] = ThetaBin;
    sgH_pre[sgH_bin] = Prelim;
    sgH_en[sgH_bin] = Energy;
    sgH_lo[sgH_bin] = EnergyLo;
    sgH_hi[sgH_bin] = EnergyHi;
    sgH_wt[sgH_bin] = Weight;
    sgH_sy[sgH_bin] = System;
    sgH_sc[sgH_bin] = Scale;
    strcpy(sgH_id[sgH_bin], Ident);

    //Increase energy bin counter
    sgH_bin++;
  }

  fclose(File_sgH);
  Sort_sgH(0, sgH_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgH_bin; t++) n+=sgH_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgH_bin; t++) if(sgH_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgH_bin);
  for(Int_t e=0; e<sgH_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgH_en[e], sgH_pts[e]);
    for(Int_t th=0; th<sgH_pts[e]; th++)
      printf("%f %f %f\n", sgH_th[e][th], sgH_val[e][th], sgH_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgH(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nH = 0;

  for(Int_t e=0; e<sgH_bin; e++)
    if((gEnergy > sgH_lo[e]) && (gEnergy < sgH_hi[e]) && (USE_PRELIMINARY || !sgH_pre[e]))
    {
      if(bins) bins[nH] = e;
      nH++;
    }

  return nH;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgH()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgH = 0.0;
  Int_t eH[EBINS];
  Int_t nH = GetEnergyBins_sgH(eH); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaH data
  for(Int_t n=0; n<nH; n++) //Process all found bins
  {
    Omega = sgH_en[eH[n]];
    for(Int_t th=0; th<sgH_pts[eH[n]]; th++) //Process all data points in current bin
    {
      Theta = sgH_th[eH[n]][th];
      Meas  = sgH_sc[eH[n]]*sgH_val[eH[n]][th]*f_obs[SIG_H];
      Error = sgH_sc[eH[n]]*sgH_err[eH[n]][th]*f_obs[SIG_H];
      Theo  = sigmaH(Theta, Omega);
      //printf("sgH: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgH+=(sgH_wt[eH[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgH;
}

//-----------------------------------------------------------------------------

void Sort_sgH(Int_t l, Int_t r) //Quicksort implementation on sgH data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgH_en[++i] < sgH_en[r]);
     while((sgH_en[--j] > sgH_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgH_lo[i],  &sgH_lo[j]);
     Swap(&sgH_en[i],  &sgH_en[j]);
     Swap(&sgH_hi[i],  &sgH_hi[j]);
     Swap(&sgH_wt[i],  &sgH_wt[j]);
     Swap(&sgH_sy[i],  &sgH_sy[j]);
     Swap(&sgH_sc[i],  &sgH_sc[j]);
     Swap(&sgH_pts[i], &sgH_pts[j]);
     Swap(&sgH_pre[i], &sgH_pre[j]);
     Swap(sgH_id[i],   sgH_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgH_val[i][n], &sgH_val[j][n]);
        Swap(&sgH_err[i][n], &sgH_err[j][n]);
        Swap(&sgH_unc[i][n], &sgH_unc[j][n]);
        Swap(&sgH_th[i][n],  &sgH_th[j][n]);
     }
   }
   Swap(&sgH_lo[i],  &sgH_lo[r]);
   Swap(&sgH_en[i],  &sgH_en[r]);
   Swap(&sgH_hi[i],  &sgH_hi[r]);
   Swap(&sgH_wt[i],  &sgH_wt[r]);
   Swap(&sgH_sy[i],  &sgH_sy[r]);
   Swap(&sgH_sc[i],  &sgH_sc[r]);
   Swap(&sgH_pts[i], &sgH_pts[r]);
   Swap(&sgH_pre[i], &sgH_pre[r]);
   Swap(sgH_id[i],   sgH_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgH_val[i][n], &sgH_val[r][n]);
      Swap(&sgH_err[i][n], &sgH_err[r][n]);
      Swap(&sgH_unc[i][n], &sgH_unc[r][n]);
      Swap(&sgH_th[i][n],  &sgH_th[r][n]);
   }

   Sort_sgH(l, i-1);
   Sort_sgH(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgH()
{
  Double_t Scale_sgH = 0.0;
  Int_t eH[EBINS];
  Int_t nH = GetEnergyBins_sgH(eH); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nH; n++) //Process all found bins
    Scale_sgH+=(f_obs[SIG_H]-1.0)*(f_obs[SIG_H]-1.0)*sgH_pts[eH[n]]/(sgH_sy[eH[n]]*sgH_sy[eH[n]]);

  return Scale_sgH;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgH()
{
  Int_t NPts_sgH = 0;
  Int_t eH[EBINS];
  Int_t nH = GetEnergyBins_sgH(eH); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nH; n++) //Process all found bins
    NPts_sgH+=sgH_pts[eH[n]];

  return NPts_sgH;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgH(FILE* File_sgH, Double_t* Theta, Double_t* sigmaH, Double_t* DsigmaH, Double_t* EsigmaH)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgH);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaH, DsigmaH, EsigmaH);
}

//-----------------------------------------------------------------------------
