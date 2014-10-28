#include "Provide_F.h"

//-----------------------------------------------------------------------------

void Load_F(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, F, DF, EF;
  Char_t Ident[256];
  FILE* File_F;

  printf("Loading   F  data from %s\n", Filename);
  File_F = fopen(Filename, "r");

  while(!feof(File_F))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_F, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_F, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_F(File_F, &Theta, &F, &DF, &EF)>=3)
    {
      F_val[F_bin][ThetaBin] = F;
      F_err[F_bin][ThetaBin] = DF;
      F_unc[F_bin][ThetaBin] = EF;
      F_th[F_bin][ThetaBin]  = Theta;
      if(DF!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    F_pts[F_bin] = ThetaBin;
    F_pre[F_bin] = Prelim;
    F_en[F_bin] = Energy;
    F_lo[F_bin] = EnergyLo;
    F_hi[F_bin] = EnergyHi;
    F_wt[F_bin] = Weight;
    F_sy[F_bin] = System;
    F_sc[F_bin] = Scale;
    strcpy(F_id[F_bin], Ident);

    //Increase energy bin counter
    F_bin++;
  }

  fclose(File_F);
  Sort_F(0, F_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<F_bin; t++) n+=F_pts[t];
  Int_t m = 0; for(Int_t t=0; t<F_bin; t++) if(F_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", F_bin);
  for(Int_t e=0; e<F_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, F_en[e], F_pts[e]);
    for(Int_t th=0; th<F_pts[e]; th++)
      printf("%f %f %f\n", F_th[e][th], F_val[e][th], F_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_F(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nF = 0;

  for(Int_t e=0; e<F_bin; e++)
    if((gEnergy > F_lo[e]) && (gEnergy < F_hi[e]) && (USE_PRELIMINARY || !F_pre[e]))
    {
      if(bins) bins[nF] = e;
      nF++;
    }

  return nF;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_F()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_F = 0.0;
  Int_t eF[EBINS];
  Int_t nF = GetEnergyBins_F(eF); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for F data
  for(Int_t n=0; n<nF; n++) //Process all found bins
  {
    Omega = F_en[eF[n]];
    for(Int_t th=0; th<F_pts[eF[n]]; th++) //Process all data points in current bin
    {
      Theta = F_th[eF[n]][th];
      Meas  = F_sc[eF[n]]*F_val[eF[n]][th]*f_obs[SIG_0];
      Error = F_sc[eF[n]]*F_err[eF[n]][th]*f_obs[SIG_0];
      Theo  = F(Theta, Omega);
      //printf("F: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_F+=(F_wt[eF[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_F;
}

//-----------------------------------------------------------------------------

void Sort_F(Int_t l, Int_t r) //Quicksort implementation on F data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(F_en[++i] < F_en[r]);
     while((F_en[--j] > F_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&F_lo[i],  &F_lo[j]);
     Swap(&F_en[i],  &F_en[j]);
     Swap(&F_hi[i],  &F_hi[j]);
     Swap(&F_wt[i],  &F_wt[j]);
     Swap(&F_sy[i],  &F_sy[j]);
     Swap(&F_pts[i], &F_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&F_val[i][n], &F_val[j][n]);
        Swap(&F_err[i][n], &F_err[j][n]);
        Swap(&F_unc[i][n], &F_unc[j][n]);
        Swap(&F_th[i][n],  &F_th[j][n]);
     }
   }
   Swap(&F_lo[i],  &F_lo[r]);
   Swap(&F_en[i],  &F_en[r]);
   Swap(&F_hi[i],  &F_hi[r]);
   Swap(&F_wt[i],  &F_wt[r]);
   Swap(&F_sy[i],  &F_sy[r]);
   Swap(&F_pts[i], &F_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&F_val[i][n], &F_val[r][n]);
      Swap(&F_err[i][n], &F_err[r][n]);
      Swap(&F_unc[i][n], &F_unc[r][n]);
      Swap(&F_th[i][n],  &F_th[r][n]);
   }

   Sort_F(l, i-1);
   Sort_F(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_F()
{
  Double_t Scale_F = 0.0;
  Int_t eF[EBINS];
  Int_t nF = GetEnergyBins_F(eF); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nF; n++) //Process all found bins
    Scale_F+=(1.0*F_pts[eF[n]])*(f_obs[SIG_0]-1.0)*(f_obs[SIG_0]-1.0)/(F_sy[eF[n]]*F_sy[eF[n]]);

  return Scale_F;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_F()
{
  Int_t NPts_F = 0;
  Int_t eF[EBINS];
  Int_t nF = GetEnergyBins_F(eF); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nF; n++) //Process all found bins
    NPts_F+=F_pts[eF[n]];

  return NPts_F;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_F(FILE* File_F, Double_t* Theta, Double_t* F, Double_t* DF, Double_t* EF)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_F);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, F, DF, EF);
}

//-----------------------------------------------------------------------------
