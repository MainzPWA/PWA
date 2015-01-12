#include "Provide_Oz.h"

//-----------------------------------------------------------------------------

void Load_Oz(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Oz, DOz, EOz;
  Char_t Ident[256];
  FILE* File_Oz;

  printf("Loading   Oz data from %s\n", Filename);
  File_Oz = fopen(Filename, "r");

  while(!feof(File_Oz))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Oz, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Oz, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Oz(File_Oz, &Theta, &Oz, &DOz, &EOz)>=3)
    {
      Oz_val[Oz_bin][ThetaBin] = Oz;
      Oz_err[Oz_bin][ThetaBin] = DOz;
      Oz_unc[Oz_bin][ThetaBin] = EOz;
      Oz_th[Oz_bin][ThetaBin]  = Theta;
      if(DOz!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Oz_pts[Oz_bin] = ThetaBin;
    Oz_pre[Oz_bin] = Prelim;
    Oz_en[Oz_bin] = Energy;
    Oz_lo[Oz_bin] = EnergyLo;
    Oz_hi[Oz_bin] = EnergyHi;
    Oz_wt[Oz_bin] = Weight;
    Oz_sy[Oz_bin] = System;
    Oz_sc[Oz_bin] = Scale;
    strcpy(Oz_id[Oz_bin], Ident);

    //Increase energy bin counter
    Oz_bin++;
  }

  fclose(File_Oz);
  Sort_Oz(0, Oz_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Oz_bin; t++) n+=Oz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Oz_bin; t++) if(Oz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Oz_bin);
  for(Int_t e=0; e<Oz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Oz_en[e], Oz_pts[e]);
    for(Int_t th=0; th<Oz_pts[e]; th++)
      printf("%f %f %f\n", Oz_th[e][th], Oz_val[e][th], Oz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_Oz(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nOz = 0;

  for(Int_t e=0; e<Oz_bin; e++)
    if((gEnergy > Oz_lo[e]) && (gEnergy < Oz_hi[e]) && (USE_PRELIMINARY || !Oz_pre[e]))
    {
      if(bins) bins[nOz] = e;
      nOz++;
    }

  return nOz;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Oz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Oz = 0.0;
  Int_t eOz[EBINS];
  Int_t nOz = GetEnergyBins_Oz(eOz); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for Oz data
  for(Int_t n=0; n<nOz; n++) //Process all found bins
  {
    Omega = Oz_en[eOz[n]];
    for(Int_t th=0; th<Oz_pts[eOz[n]]; th++) //Process all data points in current bin
    {
      Theta = Oz_th[eOz[n]][th];
      Meas  = Oz_sc[eOz[n]]*Oz_val[eOz[n]][th]*f_obs[ASY_OZ];
      Error = Oz_sc[eOz[n]]*Oz_err[eOz[n]][th]*f_obs[ASY_OZ];
      Theo  = Oz(Theta, Omega);
      //printf("Oz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_Oz+=(Oz_wt[eOz[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_Oz;
}

//-----------------------------------------------------------------------------

void Sort_Oz(Int_t l, Int_t r) //Quicksort implementation on Oz data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(Oz_en[++i] < Oz_en[r]);
     while((Oz_en[--j] > Oz_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&Oz_lo[i],  &Oz_lo[j]);
     Swap(&Oz_en[i],  &Oz_en[j]);
     Swap(&Oz_hi[i],  &Oz_hi[j]);
     Swap(&Oz_wt[i],  &Oz_wt[j]);
     Swap(&Oz_sy[i],  &Oz_sy[j]);
     Swap(&Oz_pts[i], &Oz_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&Oz_val[i][n], &Oz_val[j][n]);
        Swap(&Oz_err[i][n], &Oz_err[j][n]);
        Swap(&Oz_unc[i][n], &Oz_unc[j][n]);
        Swap(&Oz_th[i][n],  &Oz_th[j][n]);
     }
   }
   Swap(&Oz_lo[i],  &Oz_lo[r]);
   Swap(&Oz_en[i],  &Oz_en[r]);
   Swap(&Oz_hi[i],  &Oz_hi[r]);
   Swap(&Oz_wt[i],  &Oz_wt[r]);
   Swap(&Oz_sy[i],  &Oz_sy[r]);
   Swap(&Oz_pts[i], &Oz_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&Oz_val[i][n], &Oz_val[r][n]);
      Swap(&Oz_err[i][n], &Oz_err[r][n]);
      Swap(&Oz_unc[i][n], &Oz_unc[r][n]);
      Swap(&Oz_th[i][n],  &Oz_th[r][n]);
   }

   Sort_Oz(l, i-1);
   Sort_Oz(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_Oz()
{
  Double_t Scale_Oz = 0.0;
  Int_t eOz[EBINS];
  Int_t nOz = GetEnergyBins_Oz(eOz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nOz; n++) //Process all found bins
    Scale_Oz+=(1.0*Oz_pts[eOz[n]])*(f_obs[ASY_OZ]-1.0)*(f_obs[ASY_OZ]-1.0)/(Oz_sy[eOz[n]]*Oz_sy[eOz[n]]);

  return Scale_Oz;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_Oz()
{
  Int_t NPts_Oz = 0;
  Int_t eOz[EBINS];
  Int_t nOz = GetEnergyBins_Oz(eOz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nOz; n++) //Process all found bins
    NPts_Oz+=Oz_pts[eOz[n]];

  return NPts_Oz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Oz(FILE* File_Oz, Double_t* Theta, Double_t* Oz, Double_t* DOz, Double_t* EOz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Oz);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Oz, DOz, EOz);
}

//-----------------------------------------------------------------------------
