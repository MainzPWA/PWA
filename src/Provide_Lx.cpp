#include "Provide_Lx.h"

//-----------------------------------------------------------------------------

void Load_Lx(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Lx, DLx, ELx;
  Char_t Ident[256];
  FILE* File_Lx;

  printf("Loading   Lx data from %s\n", Filename);
  File_Lx = fopen(Filename, "r");

  while(!feof(File_Lx))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Lx, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Lx, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Lx(File_Lx, &Theta, &Lx, &DLx, &ELx)>=3)
    {
      Lx_val[Lx_bin][ThetaBin] = Lx;
      Lx_err[Lx_bin][ThetaBin] = DLx;
      Lx_unc[Lx_bin][ThetaBin] = ELx;
      Lx_th[Lx_bin][ThetaBin]  = Theta;
      if(DLx!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Lx_pts[Lx_bin] = ThetaBin;
    Lx_pre[Lx_bin] = Prelim;
    Lx_en[Lx_bin] = Energy;
    Lx_lo[Lx_bin] = EnergyLo;
    Lx_hi[Lx_bin] = EnergyHi;
    Lx_wt[Lx_bin] = Weight;
    Lx_sy[Lx_bin] = System;
    Lx_sc[Lx_bin] = Scale;
    strcpy(Lx_id[Lx_bin], Ident);

    //Increase energy bin counter
    Lx_bin++;
  }

  fclose(File_Lx);
  Sort_Lx(0, Lx_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Lx_bin; t++) n+=Lx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Lx_bin; t++) if(Lx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Lx_bin);
  for(Int_t e=0; e<Lx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Lx_en[e], Lx_pts[e]);
    for(Int_t th=0; th<Lx_pts[e]; th++)
      printf("%f %f %f\n", Lx_th[e][th], Lx_val[e][th], Lx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_Lx(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nLx = 0;

  for(Int_t e=0; e<Lx_bin; e++)
    if((gEnergy > Lx_lo[e]) && (gEnergy < Lx_hi[e]) && (USE_PRELIMINARY || !Lx_pre[e]))
    {
      if(bins) bins[nLx] = e;
      nLx++;
    }

  return nLx;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Lx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Lx = 0.0;
  Int_t eLx[EBINS];
  Int_t nLx = GetEnergyBins_Lx(eLx); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for Lx data
  for(Int_t n=0; n<nLx; n++) //Process all found bins
  {
    Omega = Lx_en[eLx[n]];
    for(Int_t th=0; th<Lx_pts[eLx[n]]; th++) //Process all data points in current bin
    {
      Theta = Lx_th[eLx[n]][th];
      Meas  = Lx_sc[eLx[n]]*Lx_val[eLx[n]][th]*f_obs[ASY_LX];
      Error = Lx_sc[eLx[n]]*Lx_err[eLx[n]][th]*f_obs[ASY_LX];
      Theo  = Lx(Theta, Omega);
      //printf("Lx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_Lx+=(Lx_wt[eLx[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_Lx;
}

//-----------------------------------------------------------------------------

void Sort_Lx(Int_t l, Int_t r) //Quicksort implementation on Lx data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(Lx_en[++i] < Lx_en[r]);
     while((Lx_en[--j] > Lx_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&Lx_lo[i],  &Lx_lo[j]);
     Swap(&Lx_en[i],  &Lx_en[j]);
     Swap(&Lx_hi[i],  &Lx_hi[j]);
     Swap(&Lx_wt[i],  &Lx_wt[j]);
     Swap(&Lx_sy[i],  &Lx_sy[j]);
     Swap(&Lx_pts[i], &Lx_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&Lx_val[i][n], &Lx_val[j][n]);
        Swap(&Lx_err[i][n], &Lx_err[j][n]);
        Swap(&Lx_unc[i][n], &Lx_unc[j][n]);
        Swap(&Lx_th[i][n],  &Lx_th[j][n]);
     }
   }
   Swap(&Lx_lo[i],  &Lx_lo[r]);
   Swap(&Lx_en[i],  &Lx_en[r]);
   Swap(&Lx_hi[i],  &Lx_hi[r]);
   Swap(&Lx_wt[i],  &Lx_wt[r]);
   Swap(&Lx_sy[i],  &Lx_sy[r]);
   Swap(&Lx_pts[i], &Lx_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&Lx_val[i][n], &Lx_val[r][n]);
      Swap(&Lx_err[i][n], &Lx_err[r][n]);
      Swap(&Lx_unc[i][n], &Lx_unc[r][n]);
      Swap(&Lx_th[i][n],  &Lx_th[r][n]);
   }

   Sort_Lx(l, i-1);
   Sort_Lx(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_Lx()
{
  Double_t Scale_Lx = 0.0;
  Int_t eLx[EBINS];
  Int_t nLx = GetEnergyBins_Lx(eLx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nLx; n++) //Process all found bins
    Scale_Lx+=(1.0*Lx_pts[eLx[n]])*(f_obs[ASY_LX]-1.0)*(f_obs[ASY_LX]-1.0)/(Lx_sy[eLx[n]]*Lx_sy[eLx[n]]);

  return Scale_Lx;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_Lx()
{
  Int_t NPts_Lx = 0;
  Int_t eLx[EBINS];
  Int_t nLx = GetEnergyBins_Lx(eLx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nLx; n++) //Process all found bins
    NPts_Lx+=Lx_pts[eLx[n]];

  return NPts_Lx;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Lx(FILE* File_Lx, Double_t* Theta, Double_t* Lx, Double_t* DLx, Double_t* ELx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Lx);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Lx, DLx, ELx);
}

//-----------------------------------------------------------------------------
