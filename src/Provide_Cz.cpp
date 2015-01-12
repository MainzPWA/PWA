#include "Provide_Cz.h"

//-----------------------------------------------------------------------------

void Load_Cz(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Cz, DCz, ECz;
  Char_t Ident[256];
  FILE* File_Cz;

  printf("Loading   Cz data from %s\n", Filename);
  File_Cz = fopen(Filename, "r");

  while(!feof(File_Cz))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Cz, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Cz, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Cz(File_Cz, &Theta, &Cz, &DCz, &ECz)>=3)
    {
      Cz_val[Cz_bin][ThetaBin] = Cz;
      Cz_err[Cz_bin][ThetaBin] = DCz;
      Cz_unc[Cz_bin][ThetaBin] = ECz;
      Cz_th[Cz_bin][ThetaBin]  = Theta;
      if(DCz!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Cz_pts[Cz_bin] = ThetaBin;
    Cz_pre[Cz_bin] = Prelim;
    Cz_en[Cz_bin] = Energy;
    Cz_lo[Cz_bin] = EnergyLo;
    Cz_hi[Cz_bin] = EnergyHi;
    Cz_wt[Cz_bin] = Weight;
    Cz_sy[Cz_bin] = System;
    Cz_sc[Cz_bin] = Scale;
    strcpy(Cz_id[Cz_bin], Ident);

    //Increase energy bin counter
    Cz_bin++;
  }

  fclose(File_Cz);
  Sort_Cz(0, Cz_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Cz_bin; t++) n+=Cz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Cz_bin; t++) if(Cz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Cz_bin);
  for(Int_t e=0; e<Cz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Cz_en[e], Cz_pts[e]);
    for(Int_t th=0; th<Cz_pts[e]; th++)
      printf("%f %f %f\n", Cz_th[e][th], Cz_val[e][th], Cz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_Cz(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nCz = 0;

  for(Int_t e=0; e<Cz_bin; e++)
    if((gEnergy > Cz_lo[e]) && (gEnergy < Cz_hi[e]) && (USE_PRELIMINARY || !Cz_pre[e]))
    {
      if(bins) bins[nCz] = e;
      nCz++;
    }

  return nCz;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Cz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Cz = 0.0;
  Int_t eCz[EBINS];
  Int_t nCz = GetEnergyBins_Cz(eCz); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for Cz data
  for(Int_t n=0; n<nCz; n++) //Process all found bins
  {
    Omega = Cz_en[eCz[n]];
    for(Int_t th=0; th<Cz_pts[eCz[n]]; th++) //Process all data points in current bin
    {
      Theta = Cz_th[eCz[n]][th];
      Meas  = Cz_sc[eCz[n]]*Cz_val[eCz[n]][th]*f_obs[ASY_CZ];
      Error = Cz_sc[eCz[n]]*Cz_err[eCz[n]][th]*f_obs[ASY_CZ];
      Theo  = Cz(Theta, Omega);
      //printf("Cz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_Cz+=(Cz_wt[eCz[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_Cz;
}

//-----------------------------------------------------------------------------

void Sort_Cz(Int_t l, Int_t r) //Quicksort implementation on Cz data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(Cz_en[++i] < Cz_en[r]);
     while((Cz_en[--j] > Cz_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&Cz_lo[i],  &Cz_lo[j]);
     Swap(&Cz_en[i],  &Cz_en[j]);
     Swap(&Cz_hi[i],  &Cz_hi[j]);
     Swap(&Cz_wt[i],  &Cz_wt[j]);
     Swap(&Cz_sy[i],  &Cz_sy[j]);
     Swap(&Cz_pts[i], &Cz_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&Cz_val[i][n], &Cz_val[j][n]);
        Swap(&Cz_err[i][n], &Cz_err[j][n]);
        Swap(&Cz_unc[i][n], &Cz_unc[j][n]);
        Swap(&Cz_th[i][n],  &Cz_th[j][n]);
     }
   }
   Swap(&Cz_lo[i],  &Cz_lo[r]);
   Swap(&Cz_en[i],  &Cz_en[r]);
   Swap(&Cz_hi[i],  &Cz_hi[r]);
   Swap(&Cz_wt[i],  &Cz_wt[r]);
   Swap(&Cz_sy[i],  &Cz_sy[r]);
   Swap(&Cz_pts[i], &Cz_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&Cz_val[i][n], &Cz_val[r][n]);
      Swap(&Cz_err[i][n], &Cz_err[r][n]);
      Swap(&Cz_unc[i][n], &Cz_unc[r][n]);
      Swap(&Cz_th[i][n],  &Cz_th[r][n]);
   }

   Sort_Cz(l, i-1);
   Sort_Cz(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_Cz()
{
  Double_t Scale_Cz = 0.0;
  Int_t eCz[EBINS];
  Int_t nCz = GetEnergyBins_Cz(eCz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nCz; n++) //Process all found bins
    Scale_Cz+=(1.0*Cz_pts[eCz[n]])*(f_obs[ASY_CZ]-1.0)*(f_obs[ASY_CZ]-1.0)/(Cz_sy[eCz[n]]*Cz_sy[eCz[n]]);

  return Scale_Cz;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_Cz()
{
  Int_t NPts_Cz = 0;
  Int_t eCz[EBINS];
  Int_t nCz = GetEnergyBins_Cz(eCz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nCz; n++) //Process all found bins
    NPts_Cz+=Cz_pts[eCz[n]];

  return NPts_Cz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Cz(FILE* File_Cz, Double_t* Theta, Double_t* Cz, Double_t* DCz, Double_t* ECz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Cz);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Cz, DCz, ECz);
}

//-----------------------------------------------------------------------------
