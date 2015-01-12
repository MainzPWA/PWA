#include "Provide_Ox.h"

//-----------------------------------------------------------------------------

void Load_Ox(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Ox, DOx, EOx;
  Char_t Ident[256];
  FILE* File_Ox;

  printf("Loading   Ox data from %s\n", Filename);
  File_Ox = fopen(Filename, "r");

  while(!feof(File_Ox))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Ox, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Ox, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Ox(File_Ox, &Theta, &Ox, &DOx, &EOx)>=3)
    {
      Ox_val[Ox_bin][ThetaBin] = Ox;
      Ox_err[Ox_bin][ThetaBin] = DOx;
      Ox_unc[Ox_bin][ThetaBin] = EOx;
      Ox_th[Ox_bin][ThetaBin]  = Theta;
      if(DOx!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Ox_pts[Ox_bin] = ThetaBin;
    Ox_pre[Ox_bin] = Prelim;
    Ox_en[Ox_bin] = Energy;
    Ox_lo[Ox_bin] = EnergyLo;
    Ox_hi[Ox_bin] = EnergyHi;
    Ox_wt[Ox_bin] = Weight;
    Ox_sy[Ox_bin] = System;
    Ox_sc[Ox_bin] = Scale;
    strcpy(Ox_id[Ox_bin], Ident);

    //Increase energy bin counter
    Ox_bin++;
  }

  fclose(File_Ox);
  Sort_Ox(0, Ox_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Ox_bin; t++) n+=Ox_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Ox_bin; t++) if(Ox_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Ox_bin);
  for(Int_t e=0; e<Ox_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Ox_en[e], Ox_pts[e]);
    for(Int_t th=0; th<Ox_pts[e]; th++)
      printf("%f %f %f\n", Ox_th[e][th], Ox_val[e][th], Ox_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_Ox(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nOx = 0;

  for(Int_t e=0; e<Ox_bin; e++)
    if((gEnergy > Ox_lo[e]) && (gEnergy < Ox_hi[e]) && (USE_PRELIMINARY || !Ox_pre[e]))
    {
      if(bins) bins[nOx] = e;
      nOx++;
    }

  return nOx;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Ox()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Ox = 0.0;
  Int_t eOx[EBINS];
  Int_t nOx = GetEnergyBins_Ox(eOx); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for Ox data
  for(Int_t n=0; n<nOx; n++) //Process all found bins
  {
    Omega = Ox_en[eOx[n]];
    for(Int_t th=0; th<Ox_pts[eOx[n]]; th++) //Process all data points in current bin
    {
      Theta = Ox_th[eOx[n]][th];
      Meas  = Ox_sc[eOx[n]]*Ox_val[eOx[n]][th]*f_obs[ASY_OX];
      Error = Ox_sc[eOx[n]]*Ox_err[eOx[n]][th]*f_obs[ASY_OX];
      Theo  = Ox(Theta, Omega);
      //printf("Ox: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_Ox+=(Ox_wt[eOx[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_Ox;
}

//-----------------------------------------------------------------------------

void Sort_Ox(Int_t l, Int_t r) //Quicksort implementation on Ox data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(Ox_en[++i] < Ox_en[r]);
     while((Ox_en[--j] > Ox_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&Ox_lo[i],  &Ox_lo[j]);
     Swap(&Ox_en[i],  &Ox_en[j]);
     Swap(&Ox_hi[i],  &Ox_hi[j]);
     Swap(&Ox_wt[i],  &Ox_wt[j]);
     Swap(&Ox_sy[i],  &Ox_sy[j]);
     Swap(&Ox_pts[i], &Ox_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&Ox_val[i][n], &Ox_val[j][n]);
        Swap(&Ox_err[i][n], &Ox_err[j][n]);
        Swap(&Ox_unc[i][n], &Ox_unc[j][n]);
        Swap(&Ox_th[i][n],  &Ox_th[j][n]);
     }
   }
   Swap(&Ox_lo[i],  &Ox_lo[r]);
   Swap(&Ox_en[i],  &Ox_en[r]);
   Swap(&Ox_hi[i],  &Ox_hi[r]);
   Swap(&Ox_wt[i],  &Ox_wt[r]);
   Swap(&Ox_sy[i],  &Ox_sy[r]);
   Swap(&Ox_pts[i], &Ox_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&Ox_val[i][n], &Ox_val[r][n]);
      Swap(&Ox_err[i][n], &Ox_err[r][n]);
      Swap(&Ox_unc[i][n], &Ox_unc[r][n]);
      Swap(&Ox_th[i][n],  &Ox_th[r][n]);
   }

   Sort_Ox(l, i-1);
   Sort_Ox(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_Ox()
{
  Double_t Scale_Ox = 0.0;
  Int_t eOx[EBINS];
  Int_t nOx = GetEnergyBins_Ox(eOx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nOx; n++) //Process all found bins
    Scale_Ox+=(1.0*Ox_pts[eOx[n]])*(f_obs[ASY_OX]-1.0)*(f_obs[ASY_OX]-1.0)/(Ox_sy[eOx[n]]*Ox_sy[eOx[n]]);

  return Scale_Ox;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_Ox()
{
  Int_t NPts_Ox = 0;
  Int_t eOx[EBINS];
  Int_t nOx = GetEnergyBins_Ox(eOx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nOx; n++) //Process all found bins
    NPts_Ox+=Ox_pts[eOx[n]];

  return NPts_Ox;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Ox(FILE* File_Ox, Double_t* Theta, Double_t* Ox, Double_t* DOx, Double_t* EOx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Ox);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Ox, DOx, EOx);
}

//-----------------------------------------------------------------------------
