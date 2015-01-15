#include "Provide_E.h"

//-----------------------------------------------------------------------------

void Load_E(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, E, DE, EE;
  Char_t Ident[256];
  FILE* File_E;

  printf("Loading   E  data from %s\n", Filename);
  File_E = fopen(Filename, "r");

  while(!feof(File_E))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_E, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_E, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_E(File_E, &Theta, &E, &DE, &EE)>=3)
    {
      E_val[E_bin][ThetaBin] = E;
      E_err[E_bin][ThetaBin] = DE;
      E_unc[E_bin][ThetaBin] = EE;
      E_th[E_bin][ThetaBin]  = Theta;
      if(DE!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    E_pts[E_bin] = ThetaBin;
    E_pre[E_bin] = Prelim;
    E_en[E_bin] = Energy;
    E_lo[E_bin] = EnergyLo;
    E_hi[E_bin] = EnergyHi;
    E_wt[E_bin] = Weight;
    E_sy[E_bin] = System;
    E_sc[E_bin] = Scale;
    strcpy(E_id[E_bin], Ident);

    //Increase energy bin counter
    E_bin++;
  }

  fclose(File_E);
  Sort_E(0, E_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<E_bin; t++) n+=E_pts[t];
  Int_t m = 0; for(Int_t t=0; t<E_bin; t++) if(E_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", E_bin);
  for(Int_t e=0; e<E_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, E_en[e], E_pts[e]);
    for(Int_t th=0; th<E_pts[e]; th++)
      printf("%f %f %f\n", E_th[e][th], E_val[e][th], E_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_E(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nE = 0;

  for(Int_t e=0; e<E_bin; e++)
    if((gEnergy > E_lo[e]) && (gEnergy < E_hi[e]) && (USE_PRELIMINARY || !E_pre[e]))
    {
      if(bins) bins[nE] = e;
      nE++;
    }

  return nE;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_E()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_E = 0.0;
  Int_t eE[EBINS];
  Int_t nE = GetEnergyBins_E(eE); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for E data
  for(Int_t n=0; n<nE; n++) //Process all found bins
  {
    Omega = E_en[eE[n]];
    for(Int_t th=0; th<E_pts[eE[n]]; th++) //Process all data points in current bin
    {
      Theta = E_th[eE[n]][th];
      Meas  = E_sc[eE[n]]*E_val[eE[n]][th]*f_obs[ASY_E];
      Error = E_sc[eE[n]]*E_err[eE[n]][th]*f_obs[ASY_E];
      Theo  = E(Theta, Omega);
      //printf("E: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_E+=(E_wt[eE[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_E;
}

//-----------------------------------------------------------------------------

void Sort_E(Int_t l, Int_t r) //Quicksort implementation on E data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(E_en[++i] < E_en[r]);
     while((E_en[--j] > E_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&E_lo[i],  &E_lo[j]);
     Swap(&E_en[i],  &E_en[j]);
     Swap(&E_hi[i],  &E_hi[j]);
     Swap(&E_wt[i],  &E_wt[j]);
     Swap(&E_sy[i],  &E_sy[j]);
     Swap(&E_sc[i],  &E_sc[j]);
     Swap(&E_pts[i], &E_pts[j]);
     Swap(&E_pre[i], &E_pre[j]);
     Swap(E_id[i],   E_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&E_val[i][n], &E_val[j][n]);
        Swap(&E_err[i][n], &E_err[j][n]);
        Swap(&E_unc[i][n], &E_unc[j][n]);
        Swap(&E_th[i][n],  &E_th[j][n]);
     }
   }
   Swap(&E_lo[i],  &E_lo[r]);
   Swap(&E_en[i],  &E_en[r]);
   Swap(&E_hi[i],  &E_hi[r]);
   Swap(&E_wt[i],  &E_wt[r]);
   Swap(&E_sy[i],  &E_sy[r]);
   Swap(&E_sc[i],  &E_sc[r]);
   Swap(&E_pts[i], &E_pts[r]);
   Swap(&E_pre[i], &E_pre[r]);
   Swap(E_id[i],   E_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&E_val[i][n], &E_val[r][n]);
      Swap(&E_err[i][n], &E_err[r][n]);
      Swap(&E_unc[i][n], &E_unc[r][n]);
      Swap(&E_th[i][n],  &E_th[r][n]);
   }

   Sort_E(l, i-1);
   Sort_E(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_E()
{
  Double_t Scale_E = 0.0;
  Int_t eE[EBINS];
  Int_t nE = GetEnergyBins_E(eE); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nE; n++) //Process all found bins
    Scale_E+=(f_obs[ASY_E]-1.0)*(f_obs[ASY_E]-1.0)*E_pts[eE[n]]/(E_sy[eE[n]]*E_sy[eE[n]]);

  return Scale_E;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_E()
{
  Int_t NPts_E = 0;
  Int_t eE[EBINS];
  Int_t nE = GetEnergyBins_E(eE); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nE; n++) //Process all found bins
    NPts_E+=E_pts[eE[n]];

  return NPts_E;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_E(FILE* File_E, Double_t* Theta, Double_t* E, Double_t* DE, Double_t* EE)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_E);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, E, DE, EE);
}

//-----------------------------------------------------------------------------
