#include "Provide_Tz.h"

//-----------------------------------------------------------------------------

void Load_Tz(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Tz, DTz, ETz;
  Char_t Ident[256];
  FILE* File_Tz;

  printf("Loading   Tz data from %s\n", Filename);
  File_Tz = fopen(Filename, "r");

  while(!feof(File_Tz))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Tz, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Tz, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Tz(File_Tz, &Theta, &Tz, &DTz, &ETz)>=3)
    {
      Tz_val[Tz_bin][ThetaBin] = Tz;
      Tz_err[Tz_bin][ThetaBin] = DTz;
      Tz_unc[Tz_bin][ThetaBin] = ETz;
      Tz_th[Tz_bin][ThetaBin]  = Theta;
      if(DTz!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Tz_pts[Tz_bin] = ThetaBin;
    Tz_pre[Tz_bin] = Prelim;
    Tz_en[Tz_bin] = Energy;
    Tz_lo[Tz_bin] = EnergyLo;
    Tz_hi[Tz_bin] = EnergyHi;
    Tz_wt[Tz_bin] = Weight;
    Tz_sy[Tz_bin] = System;
    Tz_sc[Tz_bin] = Scale;
    strcpy(Tz_id[Tz_bin], Ident);

    //Increase energy bin counter
    Tz_bin++;
  }

  fclose(File_Tz);
  Sort_Tz(0, Tz_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Tz_bin; t++) n+=Tz_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Tz_bin; t++) if(Tz_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Tz_bin);
  for(Int_t e=0; e<Tz_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Tz_en[e], Tz_pts[e]);
    for(Int_t th=0; th<Tz_pts[e]; th++)
      printf("%f %f %f\n", Tz_th[e][th], Tz_val[e][th], Tz_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_Tz(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nTz = 0;

  for(Int_t e=0; e<Tz_bin; e++)
    if((gEnergy > Tz_lo[e]) && (gEnergy < Tz_hi[e]) && (USE_PRELIMINARY || !Tz_pre[e]))
    {
      if(bins) bins[nTz] = e;
      nTz++;
    }

  return nTz;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Tz()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Tz = 0.0;
  Int_t eTz[EBINS];
  Int_t nTz = GetEnergyBins_Tz(eTz); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for Tz data
  for(Int_t n=0; n<nTz; n++) //Process all found bins
  {
    Omega = Tz_en[eTz[n]];
    for(Int_t th=0; th<Tz_pts[eTz[n]]; th++) //Process all data points in current bin
    {
      Theta = Tz_th[eTz[n]][th];
      Meas  = Tz_sc[eTz[n]]*Tz_val[eTz[n]][th]*f_obs[ASY_TZ];
      Error = Tz_sc[eTz[n]]*Tz_err[eTz[n]][th]*f_obs[ASY_TZ];
      Theo  = Tz(Theta, Omega);
      //printf("Tz: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_Tz+=(Tz_wt[eTz[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_Tz;
}

//-----------------------------------------------------------------------------

void Sort_Tz(Int_t l, Int_t r) //Quicksort implementation on Tz data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(Tz_en[++i] < Tz_en[r]);
     while((Tz_en[--j] > Tz_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&Tz_lo[i],  &Tz_lo[j]);
     Swap(&Tz_en[i],  &Tz_en[j]);
     Swap(&Tz_hi[i],  &Tz_hi[j]);
     Swap(&Tz_wt[i],  &Tz_wt[j]);
     Swap(&Tz_sy[i],  &Tz_sy[j]);
     Swap(&Tz_sc[i],  &Tz_sc[j]);
     Swap(&Tz_pts[i], &Tz_pts[j]);
     Swap(&Tz_pre[i], &Tz_pre[j]);
     Swap(Tz_id[i],   Tz_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&Tz_val[i][n], &Tz_val[j][n]);
        Swap(&Tz_err[i][n], &Tz_err[j][n]);
        Swap(&Tz_unc[i][n], &Tz_unc[j][n]);
        Swap(&Tz_th[i][n],  &Tz_th[j][n]);
     }
   }
   Swap(&Tz_lo[i],  &Tz_lo[r]);
   Swap(&Tz_en[i],  &Tz_en[r]);
   Swap(&Tz_hi[i],  &Tz_hi[r]);
   Swap(&Tz_wt[i],  &Tz_wt[r]);
   Swap(&Tz_sy[i],  &Tz_sy[r]);
   Swap(&Tz_sc[i],  &Tz_sc[r]);
   Swap(&Tz_pts[i], &Tz_pts[r]);
   Swap(&Tz_pre[i], &Tz_pre[r]);
   Swap(Tz_id[i],   Tz_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&Tz_val[i][n], &Tz_val[r][n]);
      Swap(&Tz_err[i][n], &Tz_err[r][n]);
      Swap(&Tz_unc[i][n], &Tz_unc[r][n]);
      Swap(&Tz_th[i][n],  &Tz_th[r][n]);
   }

   Sort_Tz(l, i-1);
   Sort_Tz(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_Tz()
{
  Double_t Scale_Tz = 0.0;
  Int_t eTz[EBINS];
  Int_t nTz = GetEnergyBins_Tz(eTz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nTz; n++) //Process all found bins
    Scale_Tz+=(f_obs[ASY_TZ]-1.0)*(f_obs[ASY_TZ]-1.0)*Tz_pts[eTz[n]]/(Tz_sy[eTz[n]]*Tz_sy[eTz[n]]);

  return Scale_Tz;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_Tz()
{
  Int_t NPts_Tz = 0;
  Int_t eTz[EBINS];
  Int_t nTz = GetEnergyBins_Tz(eTz); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nTz; n++) //Process all found bins
    NPts_Tz+=Tz_pts[eTz[n]];

  return NPts_Tz;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Tz(FILE* File_Tz, Double_t* Theta, Double_t* Tz, Double_t* DTz, Double_t* ETz)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Tz);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Tz, DTz, ETz);
}

//-----------------------------------------------------------------------------
