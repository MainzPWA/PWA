#include "Provide_H.h"

//-----------------------------------------------------------------------------

void Load_H(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, H, DH, EH;
  Char_t Ident[256];
  FILE* File_H;

  printf("Loading   H  data from %s\n", Filename);
  File_H = fopen(Filename, "r");

  while(!feof(File_H))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_H, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_H, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_H(File_H, &Theta, &H, &DH, &EH)>=3)
    {
      H_val[H_bin][ThetaBin] = H;
      H_err[H_bin][ThetaBin] = DH;
      H_unc[H_bin][ThetaBin] = EH;
      H_th[H_bin][ThetaBin]  = Theta;
      if(DH!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    H_pts[H_bin] = ThetaBin;
    H_pre[H_bin] = Prelim;
    H_en[H_bin] = Energy;
    H_lo[H_bin] = EnergyLo;
    H_hi[H_bin] = EnergyHi;
    H_wt[H_bin] = Weight;
    H_sy[H_bin] = System;
    H_sc[H_bin] = Scale;
    strcpy(H_id[H_bin], Ident);

    //Increase energy bin counter
    H_bin++;
  }

  fclose(File_H);
  Sort_H(0, H_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<H_bin; t++) n+=H_pts[t];
  Int_t m = 0; for(Int_t t=0; t<H_bin; t++) if(H_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", H_bin);
  for(Int_t e=0; e<H_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, H_en[e], H_pts[e]);
    for(Int_t th=0; th<H_pts[e]; th++)
      printf("%f %f %f\n", H_th[e][th], H_val[e][th], H_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_H(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nH = 0;

  for(Int_t e=0; e<H_bin; e++)
    if((gEnergy > H_lo[e]) && (gEnergy < H_hi[e]) && (USE_PRELIMINARY || !H_pre[e]))
    {
      if(bins) bins[nH] = e;
      nH++;
    }

  return nH;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_H()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_H = 0.0;
  Int_t eH[EBINS];
  Int_t nH = GetEnergyBins_H(eH); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for H data
  for(Int_t n=0; n<nH; n++) //Process all found bins
  {
    Omega = H_en[eH[n]];
    for(Int_t th=0; th<H_pts[eH[n]]; th++) //Process all data points in current bin
    {
      Theta = H_th[eH[n]][th];
      Meas  = H_sc[eH[n]]*H_val[eH[n]][th]*f_obs[SIG_0];
      Error = H_sc[eH[n]]*H_err[eH[n]][th]*f_obs[SIG_0];
      Theo  = H(Theta, Omega);
      //printf("H: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_H+=(H_wt[eH[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_H;
}

//-----------------------------------------------------------------------------

void Sort_H(Int_t l, Int_t r) //Quicksort implementation on H data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(H_en[++i] < H_en[r]);
     while((H_en[--j] > H_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&H_lo[i],  &H_lo[j]);
     Swap(&H_en[i],  &H_en[j]);
     Swap(&H_hi[i],  &H_hi[j]);
     Swap(&H_wt[i],  &H_wt[j]);
     Swap(&H_sy[i],  &H_sy[j]);
     Swap(&H_pts[i], &H_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&H_val[i][n], &H_val[j][n]);
        Swap(&H_err[i][n], &H_err[j][n]);
        Swap(&H_unc[i][n], &H_unc[j][n]);
        Swap(&H_th[i][n],  &H_th[j][n]);
     }
   }
   Swap(&H_lo[i],  &H_lo[r]);
   Swap(&H_en[i],  &H_en[r]);
   Swap(&H_hi[i],  &H_hi[r]);
   Swap(&H_wt[i],  &H_wt[r]);
   Swap(&H_sy[i],  &H_sy[r]);
   Swap(&H_pts[i], &H_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&H_val[i][n], &H_val[r][n]);
      Swap(&H_err[i][n], &H_err[r][n]);
      Swap(&H_unc[i][n], &H_unc[r][n]);
      Swap(&H_th[i][n],  &H_th[r][n]);
   }

   Sort_H(l, i-1);
   Sort_H(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_H()
{
  Double_t Scale_H = 0.0;
  Int_t eH[EBINS];
  Int_t nH = GetEnergyBins_H(eH); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nH; n++) //Process all found bins
    Scale_H+=(1.0*H_pts[eH[n]])*(f_obs[SIG_0]-1.0)*(f_obs[SIG_0]-1.0)/(H_sy[eH[n]]*H_sy[eH[n]]);

  return Scale_H;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_H()
{
  Int_t NPts_H = 0;
  Int_t eH[EBINS];
  Int_t nH = GetEnergyBins_H(eH); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nH; n++) //Process all found bins
    NPts_H+=H_pts[eH[n]];

  return NPts_H;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_H(FILE* File_H, Double_t* Theta, Double_t* H, Double_t* DH, Double_t* EH)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_H);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, H, DH, EH);
}

//-----------------------------------------------------------------------------
