#include "Provide_T.h"

//-----------------------------------------------------------------------------

void Load_T(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, T, DT, ET;
  Char_t Ident[256];
  FILE* File_T;

  printf("Loading   T  data from %s\n", Filename);
  File_T = fopen(Filename, "r");

  while(!feof(File_T))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_T, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_T, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_T(File_T, &Theta, &T, &DT, &ET)>=3)
    {
      T_val[T_bin][ThetaBin] = T;
      T_err[T_bin][ThetaBin] = DT;
      T_unc[T_bin][ThetaBin] = ET;
      T_th[T_bin][ThetaBin]  = Theta;
      if(DT!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    T_pts[T_bin] = ThetaBin;
    T_pre[T_bin] = Prelim;
    T_en[T_bin] = Energy;
    T_lo[T_bin] = EnergyLo;
    T_hi[T_bin] = EnergyHi;
    T_wt[T_bin] = Weight;
    T_sy[T_bin] = System;
    T_sc[T_bin] = Scale;
    strcpy(T_id[T_bin], Ident);

    //Increase energy bin counter
    T_bin++;
  }

  fclose(File_T);
  Sort_T(0, T_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<T_bin; t++) n+=T_pts[t];
  Int_t m = 0; for(Int_t t=0; t<T_bin; t++) if(T_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", T_bin);
  for(Int_t e=0; e<T_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, T_en[e], T_pts[e]);
    for(Int_t th=0; th<T_pts[e]; th++)
      printf("%f %f %f\n", T_th[e][th], T_val[e][th], T_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_T(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nT = 0;

  for(Int_t e=0; e<T_bin; e++)
    if((gEnergy > T_lo[e]) && (gEnergy < T_hi[e]) && (USE_PRELIMINARY || !T_pre[e]))
    {
      if(bins) bins[nT] = e;
      nT++;
    }

  return nT;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_T()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_T = 0.0;
  Int_t eT[EBINS];
  Int_t nT = GetEnergyBins_T(eT); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for T data
  for(Int_t n=0; n<nT; n++) //Process all found bins
  {
    Omega = T_en[eT[n]];
    for(Int_t th=0; th<T_pts[eT[n]]; th++) //Process all data points in current bin
    {
      Theta = T_th[eT[n]][th];
      Meas  = T_sc[eT[n]]*T_val[eT[n]][th]*f_obs[SIG_0];
      Error = T_sc[eT[n]]*T_err[eT[n]][th]*f_obs[SIG_0];
      Theo  = T(Theta, Omega);
      printf("T: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_T+=(T_wt[eT[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_T;
}

//-----------------------------------------------------------------------------

void Sort_T(Int_t l, Int_t r) //Quicksort implementation on T data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(T_en[++i] < T_en[r]);
     while((T_en[--j] > T_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&T_lo[i],  &T_lo[j]);
     Swap(&T_en[i],  &T_en[j]);
     Swap(&T_hi[i],  &T_hi[j]);
     Swap(&T_wt[i],  &T_wt[j]);
     Swap(&T_sy[i],  &T_sy[j]);
     Swap(&T_pts[i], &T_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&T_val[i][n], &T_val[j][n]);
        Swap(&T_err[i][n], &T_err[j][n]);
        Swap(&T_unc[i][n], &T_unc[j][n]);
        Swap(&T_th[i][n],  &T_th[j][n]);
     }
   }
   Swap(&T_lo[i],  &T_lo[r]);
   Swap(&T_en[i],  &T_en[r]);
   Swap(&T_hi[i],  &T_hi[r]);
   Swap(&T_wt[i],  &T_wt[r]);
   Swap(&T_sy[i],  &T_sy[r]);
   Swap(&T_pts[i], &T_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&T_val[i][n], &T_val[r][n]);
      Swap(&T_err[i][n], &T_err[r][n]);
      Swap(&T_unc[i][n], &T_unc[r][n]);
      Swap(&T_th[i][n],  &T_th[r][n]);
   }

   Sort_T(l, i-1);
   Sort_T(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_T()
{
  Double_t Scale_T = 0.0;
  Int_t eT[EBINS];
  Int_t nT = GetEnergyBins_T(eT); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nT; n++) //Process all found bins
    Scale_T+=(1.0*T_pts[eT[n]])*(f_obs[SIG_0]-1.0)*(f_obs[SIG_0]-1.0)/(T_sy[eT[n]]*T_sy[eT[n]]);

  return Scale_T;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_T()
{
  Int_t NPts_T = 0;
  Int_t eT[EBINS];
  Int_t nT = GetEnergyBins_T(eT); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nT; n++) //Process all found bins
    NPts_T+=T_pts[eT[n]];

  return NPts_T;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_T(FILE* File_T, Double_t* Theta, Double_t* T, Double_t* DT, Double_t* ET)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_T);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, T, DT, ET);
}

//-----------------------------------------------------------------------------
