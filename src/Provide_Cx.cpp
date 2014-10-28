#include "Provide_Cx.h"

//-----------------------------------------------------------------------------

void Load_Cx(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Cx, DCx, ECx;
  Char_t Ident[256];
  FILE* File_Cx;

  printf("Loading   Cx data from %s\n", Filename);
  File_Cx = fopen(Filename, "r");

  while(!feof(File_Cx))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Cx, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Cx, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Cx(File_Cx, &Theta, &Cx, &DCx, &ECx)>=3)
    {
      Cx_val[Cx_bin][ThetaBin] = Cx;
      Cx_err[Cx_bin][ThetaBin] = DCx;
      Cx_unc[Cx_bin][ThetaBin] = ECx;
      Cx_th[Cx_bin][ThetaBin]  = Theta;
      if(DCx!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Cx_pts[Cx_bin] = ThetaBin;
    Cx_pre[Cx_bin] = Prelim;
    Cx_en[Cx_bin] = Energy;
    Cx_lo[Cx_bin] = EnergyLo;
    Cx_hi[Cx_bin] = EnergyHi;
    Cx_wt[Cx_bin] = Weight;
    Cx_sy[Cx_bin] = System;
    Cx_sc[Cx_bin] = Scale;
    strcpy(Cx_id[Cx_bin], Ident);

    //Increase energy bin counter
    Cx_bin++;
  }

  fclose(File_Cx);
  Sort_Cx(0, Cx_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Cx_bin; t++) n+=Cx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Cx_bin; t++) if(Cx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Cx_bin);
  for(Int_t e=0; e<Cx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Cx_en[e], Cx_pts[e]);
    for(Int_t th=0; th<Cx_pts[e]; th++)
      printf("%f %f %f\n", Cx_th[e][th], Cx_val[e][th], Cx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_Cx(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nCx = 0;

  for(Int_t e=0; e<Cx_bin; e++)
    if((gEnergy > Cx_lo[e]) && (gEnergy < Cx_hi[e]) && (USE_PRELIMINARY || !Cx_pre[e]))
    {
      if(bins) bins[nCx] = e;
      nCx++;
    }

  return nCx;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Cx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Cx = 0.0;
  Int_t eCx[EBINS];
  Int_t nCx = GetEnergyBins_Cx(eCx); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for Cx data
  for(Int_t n=0; n<nCx; n++) //Process all found bins
  {
    Omega = Cx_en[eCx[n]];
    for(Int_t th=0; th<Cx_pts[eCx[n]]; th++) //Process all data points in current bin
    {
      Theta = Cx_th[eCx[n]][th];
      Meas  = Cx_sc[eCx[n]]*Cx_val[eCx[n]][th]*f_obs[SIG_0];
      Error = Cx_sc[eCx[n]]*Cx_err[eCx[n]][th]*f_obs[SIG_0];
      Theo  = Cx(Theta, Omega);
      //printf("Cx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_Cx+=(Cx_wt[eCx[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_Cx;
}

//-----------------------------------------------------------------------------

void Sort_Cx(Int_t l, Int_t r) //Quicksort implementation on Cx data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(Cx_en[++i] < Cx_en[r]);
     while((Cx_en[--j] > Cx_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&Cx_lo[i],  &Cx_lo[j]);
     Swap(&Cx_en[i],  &Cx_en[j]);
     Swap(&Cx_hi[i],  &Cx_hi[j]);
     Swap(&Cx_wt[i],  &Cx_wt[j]);
     Swap(&Cx_sy[i],  &Cx_sy[j]);
     Swap(&Cx_pts[i], &Cx_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&Cx_val[i][n], &Cx_val[j][n]);
        Swap(&Cx_err[i][n], &Cx_err[j][n]);
        Swap(&Cx_unc[i][n], &Cx_unc[j][n]);
        Swap(&Cx_th[i][n],  &Cx_th[j][n]);
     }
   }
   Swap(&Cx_lo[i],  &Cx_lo[r]);
   Swap(&Cx_en[i],  &Cx_en[r]);
   Swap(&Cx_hi[i],  &Cx_hi[r]);
   Swap(&Cx_wt[i],  &Cx_wt[r]);
   Swap(&Cx_sy[i],  &Cx_sy[r]);
   Swap(&Cx_pts[i], &Cx_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&Cx_val[i][n], &Cx_val[r][n]);
      Swap(&Cx_err[i][n], &Cx_err[r][n]);
      Swap(&Cx_unc[i][n], &Cx_unc[r][n]);
      Swap(&Cx_th[i][n],  &Cx_th[r][n]);
   }

   Sort_Cx(l, i-1);
   Sort_Cx(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_Cx()
{
  Double_t Scale_Cx = 0.0;
  Int_t eCx[EBINS];
  Int_t nCx = GetEnergyBins_Cx(eCx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nCx; n++) //Process all found bins
    Scale_Cx+=(1.0*Cx_pts[eCx[n]])*(f_obs[SIG_0]-1.0)*(f_obs[SIG_0]-1.0)/(Cx_sy[eCx[n]]*Cx_sy[eCx[n]]);

  return Scale_Cx;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_Cx()
{
  Int_t NPts_Cx = 0;
  Int_t eCx[EBINS];
  Int_t nCx = GetEnergyBins_Cx(eCx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nCx; n++) //Process all found bins
    NPts_Cx+=Cx_pts[eCx[n]];

  return NPts_Cx;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Cx(FILE* File_Cx, Double_t* Theta, Double_t* Cx, Double_t* DCx, Double_t* ECx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Cx);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Cx, DCx, ECx);
}

//-----------------------------------------------------------------------------
