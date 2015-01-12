#include "Provide_P.h"

//-----------------------------------------------------------------------------

void Load_P(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, P, DP, EP;
  Char_t Ident[256];
  FILE* File_P;

  printf("Loading   P  data from %s\n", Filename);
  File_P = fopen(Filename, "r");

  while(!feof(File_P))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_P, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_P, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_P(File_P, &Theta, &P, &DP, &EP)>=3)
    {
      P_val[P_bin][ThetaBin] = P;
      P_err[P_bin][ThetaBin] = DP;
      P_unc[P_bin][ThetaBin] = EP;
      P_th[P_bin][ThetaBin]  = Theta;
      if(DP!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    P_pts[P_bin] = ThetaBin;
    P_pre[P_bin] = Prelim;
    P_en[P_bin] = Energy;
    P_lo[P_bin] = EnergyLo;
    P_hi[P_bin] = EnergyHi;
    P_wt[P_bin] = Weight;
    P_sy[P_bin] = System;
    P_sc[P_bin] = Scale;
    strcpy(P_id[P_bin], Ident);

    //Increase energy bin counter
    P_bin++;
  }

  fclose(File_P);
  Sort_P(0, P_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<P_bin; t++) n+=P_pts[t];
  Int_t m = 0; for(Int_t t=0; t<P_bin; t++) if(P_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", P_bin);
  for(Int_t e=0; e<P_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, P_en[e], P_pts[e]);
    for(Int_t th=0; th<P_pts[e]; th++)
      printf("%f %f %f\n", P_th[e][th], P_val[e][th], P_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_P(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nP = 0;

  for(Int_t e=0; e<P_bin; e++)
    if((gEnergy > P_lo[e]) && (gEnergy < P_hi[e]) && (USE_PRELIMINARY || !P_pre[e]))
    {
      if(bins) bins[nP] = e;
      nP++;
    }

  return nP;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_P()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_P = 0.0;
  Int_t eP[EBINS];
  Int_t nP = GetEnergyBins_P(eP); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for P data
  for(Int_t n=0; n<nP; n++) //Process all found bins
  {
    Omega = P_en[eP[n]];
    for(Int_t th=0; th<P_pts[eP[n]]; th++) //Process all data points in current bin
    {
      Theta = P_th[eP[n]][th];
      Meas  = P_sc[eP[n]]*P_val[eP[n]][th]*f_obs[ASY_P];
      Error = P_sc[eP[n]]*P_err[eP[n]][th]*f_obs[ASY_P];
      Theo  = P(Theta, Omega);
      //printf("P: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_P+=(P_wt[eP[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_P;
}

//-----------------------------------------------------------------------------

void Sort_P(Int_t l, Int_t r) //Quicksort implementation on P data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(P_en[++i] < P_en[r]);
     while((P_en[--j] > P_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&P_lo[i],  &P_lo[j]);
     Swap(&P_en[i],  &P_en[j]);
     Swap(&P_hi[i],  &P_hi[j]);
     Swap(&P_wt[i],  &P_wt[j]);
     Swap(&P_sy[i],  &P_sy[j]);
     Swap(&P_pts[i], &P_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&P_val[i][n], &P_val[j][n]);
        Swap(&P_err[i][n], &P_err[j][n]);
        Swap(&P_unc[i][n], &P_unc[j][n]);
        Swap(&P_th[i][n],  &P_th[j][n]);
     }
   }
   Swap(&P_lo[i],  &P_lo[r]);
   Swap(&P_en[i],  &P_en[r]);
   Swap(&P_hi[i],  &P_hi[r]);
   Swap(&P_wt[i],  &P_wt[r]);
   Swap(&P_sy[i],  &P_sy[r]);
   Swap(&P_pts[i], &P_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&P_val[i][n], &P_val[r][n]);
      Swap(&P_err[i][n], &P_err[r][n]);
      Swap(&P_unc[i][n], &P_unc[r][n]);
      Swap(&P_th[i][n],  &P_th[r][n]);
   }

   Sort_P(l, i-1);
   Sort_P(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_P()
{
  Double_t Scale_P = 0.0;
  Int_t eP[EBINS];
  Int_t nP = GetEnergyBins_P(eP); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nP; n++) //Process all found bins
    Scale_P+=(1.0*P_pts[eP[n]])*(f_obs[ASY_P]-1.0)*(f_obs[ASY_P]-1.0)/(P_sy[eP[n]]*P_sy[eP[n]]);

  return Scale_P;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_P()
{
  Int_t NPts_P = 0;
  Int_t eP[EBINS];
  Int_t nP = GetEnergyBins_P(eP); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nP; n++) //Process all found bins
    NPts_P+=P_pts[eP[n]];

  return NPts_P;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_P(FILE* File_P, Double_t* Theta, Double_t* P, Double_t* DP, Double_t* EP)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_P);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, P, DP, EP);
}

//-----------------------------------------------------------------------------
