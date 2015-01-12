#include "Provide_Lz.h"

//-----------------------------------------------------------------------------

void Load_Lz(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Lz, DLz, ELz;
  Char_t Ident[256];
  FILE* File_Lz;

  printf("Loading   Lz data from %s\n", Filename);
  File_Lz = fopen(Filename, "r");

  while(!feof(File_Lz))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Lz, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Lz, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Lz(File_Lz, &Theta, &Lz, &DLz, &ELz)>=3)
    {
      Lz_val[Lz_bin][ThetaBin] = Lz;
      Lz_err[Lz_bin][ThetaBin] = DLz;
      Lz_unc[Lz_bin][ThetaBin] = ELz;
      Lz_th[Lz_bin][ThetaBin]  = Theta;
      if(DLz!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Lz_pts[Lz_bin] = ThetaBin;
    Lz_pre[Lz_bin] = Prelim;
    Lz_en[Lz_bin] = Energy;
    Lz_lo[Lz_bin] = EnergyLo;
    Lz_hi[Lz_bin] = EnergyHi;
    Lz_wt[Lz_bin] = Weight;
    Lz_sy[Lz_bin] = System;
    Lz_sc[Lz_bin] = Scale;
    strcpy(Lz_id[Lz_bin], Ident);

    //Increase energy bin counter
    Lz_bin++;
  }

  fclose(File_Lz);
  Sort_Lz(0, Lz_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Lz_bin; t++) n+=Lz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Lz_bin; t++) if(Lz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Lz_bin);
  for(Int_t e=0; e<Lz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Lz_en[e], Lz_pts[e]);
    for(Int_t th=0; th<Lz_pts[e]; th++)
      printf("%f %f %f\n", Lz_th[e][th], Lz_val[e][th], Lz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_Lz(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nLz = 0;

  for(Int_t e=0; e<Lz_bin; e++)
    if((gEnergy > Lz_lo[e]) && (gEnergy < Lz_hi[e]) && (USE_PRELIMINARY || !Lz_pre[e]))
    {
      if(bins) bins[nLz] = e;
      nLz++;
    }

  return nLz;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Lz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Lz = 0.0;
  Int_t eLz[EBINS];
  Int_t nLz = GetEnergyBins_Lz(eLz); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for Lz data
  for(Int_t n=0; n<nLz; n++) //Process all found bins
  {
    Omega = Lz_en[eLz[n]];
    for(Int_t th=0; th<Lz_pts[eLz[n]]; th++) //Process all data points in current bin
    {
      Theta = Lz_th[eLz[n]][th];
      Meas  = Lz_sc[eLz[n]]*Lz_val[eLz[n]][th]*f_obs[ASY_LZ];
      Error = Lz_sc[eLz[n]]*Lz_err[eLz[n]][th]*f_obs[ASY_LZ];
      Theo  = Lz(Theta, Omega);
      //printf("Lz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_Lz+=(Lz_wt[eLz[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_Lz;
}

//-----------------------------------------------------------------------------

void Sort_Lz(Int_t l, Int_t r) //Quicksort implementation on Lz data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(Lz_en[++i] < Lz_en[r]);
     while((Lz_en[--j] > Lz_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&Lz_lo[i],  &Lz_lo[j]);
     Swap(&Lz_en[i],  &Lz_en[j]);
     Swap(&Lz_hi[i],  &Lz_hi[j]);
     Swap(&Lz_wt[i],  &Lz_wt[j]);
     Swap(&Lz_sy[i],  &Lz_sy[j]);
     Swap(&Lz_pts[i], &Lz_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&Lz_val[i][n], &Lz_val[j][n]);
        Swap(&Lz_err[i][n], &Lz_err[j][n]);
        Swap(&Lz_unc[i][n], &Lz_unc[j][n]);
        Swap(&Lz_th[i][n],  &Lz_th[j][n]);
     }
   }
   Swap(&Lz_lo[i],  &Lz_lo[r]);
   Swap(&Lz_en[i],  &Lz_en[r]);
   Swap(&Lz_hi[i],  &Lz_hi[r]);
   Swap(&Lz_wt[i],  &Lz_wt[r]);
   Swap(&Lz_sy[i],  &Lz_sy[r]);
   Swap(&Lz_pts[i], &Lz_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&Lz_val[i][n], &Lz_val[r][n]);
      Swap(&Lz_err[i][n], &Lz_err[r][n]);
      Swap(&Lz_unc[i][n], &Lz_unc[r][n]);
      Swap(&Lz_th[i][n],  &Lz_th[r][n]);
   }

   Sort_Lz(l, i-1);
   Sort_Lz(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_Lz()
{
  Double_t Scale_Lz = 0.0;
  Int_t eLz[EBINS];
  Int_t nLz = GetEnergyBins_Lz(eLz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nLz; n++) //Process all found bins
    Scale_Lz+=(f_obs[ASY_LZ]-1.0)*(f_obs[ASY_LZ]-1.0)*Lz_pts[eLz[n]]/(Lz_sy[eLz[n]]*Lz_sy[eLz[n]]);

  return Scale_Lz;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_Lz()
{
  Int_t NPts_Lz = 0;
  Int_t eLz[EBINS];
  Int_t nLz = GetEnergyBins_Lz(eLz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nLz; n++) //Process all found bins
    NPts_Lz+=Lz_pts[eLz[n]];

  return NPts_Lz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Lz(FILE* File_Lz, Double_t* Theta, Double_t* Lz, Double_t* DLz, Double_t* ELz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Lz);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Lz, DLz, ELz);
}

//-----------------------------------------------------------------------------
