#include "Provide_S.h"

//-----------------------------------------------------------------------------

void Load_S(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, S, DS, ES;
  Char_t Ident[256];
  FILE* File_S;

  printf("Loading   S  data from %s\n", Filename);
  File_S = fopen(Filename, "r");

  while(!feof(File_S))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_S, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_S, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_S(File_S, &Theta, &S, &DS, &ES)>=3)
    {
      S_val[S_bin][ThetaBin] = S;
      S_err[S_bin][ThetaBin] = DS;
      S_unc[S_bin][ThetaBin] = ES;
      S_th[S_bin][ThetaBin]  = Theta;
      if(DS!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    S_pts[S_bin] = ThetaBin;
    S_pre[S_bin] = Prelim;
    S_en[S_bin] = Energy;
    S_lo[S_bin] = EnergyLo;
    S_hi[S_bin] = EnergyHi;
    S_wt[S_bin] = Weight;
    S_sy[S_bin] = System;
    S_sc[S_bin] = Scale;
    strcpy(S_id[S_bin], Ident);

    //Increase energy bin counter
    S_bin++;
  }

  fclose(File_S);
  Sort_S(0, S_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<S_bin; t++) n+=S_pts[t];
  Int_t m = 0; for(Int_t t=0; t<S_bin; t++) if(S_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", S_bin);
  for(Int_t e=0; e<S_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, S_en[e], S_pts[e]);
    for(Int_t th=0; th<S_pts[e]; th++)
      printf("%f %f %f\n", S_th[e][th], S_val[e][th], S_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_S(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nS = 0;

  for(Int_t e=0; e<S_bin; e++)
    if((gEnergy > S_lo[e]) && (gEnergy < S_hi[e]) && (USE_PRELIMINARY || !S_pre[e]))
    {
      if(bins) bins[nS] = e;
      nS++;
    }

  return nS;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_S()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_S = 0.0;
  Int_t eS[EBINS];
  Int_t nS = GetEnergyBins_S(eS); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for S data
  for(Int_t n=0; n<nS; n++) //Process all found bins
  {
    Omega = S_en[eS[n]];
    for(Int_t th=0; th<S_pts[eS[n]]; th++) //Process all data points in current bin
    {
      Theta = S_th[eS[n]][th];
      Meas  = S_sc[eS[n]]*S_val[eS[n]][th]*f_obs[ASY_S];
      Error = S_sc[eS[n]]*S_err[eS[n]][th]*f_obs[ASY_S];
      Theo  = S(Theta, Omega);
      //printf("S: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_S+=(S_wt[eS[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_S;
}

//-----------------------------------------------------------------------------

void Sort_S(Int_t l, Int_t r) //Quicksort implementation on S data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(S_en[++i] < S_en[r]);
     while((S_en[--j] > S_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&S_lo[i],  &S_lo[j]);
     Swap(&S_en[i],  &S_en[j]);
     Swap(&S_hi[i],  &S_hi[j]);
     Swap(&S_wt[i],  &S_wt[j]);
     Swap(&S_sy[i],  &S_sy[j]);
     Swap(&S_sc[i],  &S_sc[j]);
     Swap(&S_pts[i], &S_pts[j]);
     Swap(&S_pre[i], &S_pre[j]);
     Swap(S_id[i],   S_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&S_val[i][n], &S_val[j][n]);
        Swap(&S_err[i][n], &S_err[j][n]);
        Swap(&S_unc[i][n], &S_unc[j][n]);
        Swap(&S_th[i][n],  &S_th[j][n]);
     }
   }
   Swap(&S_lo[i],  &S_lo[r]);
   Swap(&S_en[i],  &S_en[r]);
   Swap(&S_hi[i],  &S_hi[r]);
   Swap(&S_wt[i],  &S_wt[r]);
   Swap(&S_sy[i],  &S_sy[r]);
   Swap(&S_sc[i],  &S_sc[r]);
   Swap(&S_pts[i], &S_pts[r]);
   Swap(&S_pre[i], &S_pre[r]);
   Swap(S_id[i],   S_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&S_val[i][n], &S_val[r][n]);
      Swap(&S_err[i][n], &S_err[r][n]);
      Swap(&S_unc[i][n], &S_unc[r][n]);
      Swap(&S_th[i][n],  &S_th[r][n]);
   }

   Sort_S(l, i-1);
   Sort_S(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_S()
{
  Double_t Scale_S = 0.0;
  Int_t eS[EBINS];
  Int_t nS = GetEnergyBins_S(eS); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nS; n++) //Process all found bins
    Scale_S+=(f_obs[ASY_S]-1.0)*(f_obs[ASY_S]-1.0)*S_pts[eS[n]]/(S_sy[eS[n]]*S_sy[eS[n]]);

  return Scale_S;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_S()
{
  Int_t NPts_S = 0;
  Int_t eS[EBINS];
  Int_t nS = GetEnergyBins_S(eS); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nS; n++) //Process all found bins
    NPts_S+=S_pts[eS[n]];

  return NPts_S;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_S(FILE* File_S, Double_t* Theta, Double_t* S, Double_t* DS, Double_t* ES)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_S);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, S, DS, ES);
}

//-----------------------------------------------------------------------------
