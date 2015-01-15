#include "Provide_sg0.h"

//-----------------------------------------------------------------------------

void Load_sg0(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigma0, Dsigma0, Esigma0;
  Char_t Ident[256];
  FILE* File_sg0;

  printf("Loading sg0  data from %s\n", Filename);
  File_sg0 = fopen(Filename, "r");

  while(!feof(File_sg0))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sg0, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sg0, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sg0(File_sg0, &Theta, &sigma0, &Dsigma0, &Esigma0)>=3)
    {
      sg0_val[sg0_bin][ThetaBin] = sigma0;
      sg0_err[sg0_bin][ThetaBin] = Dsigma0;
      sg0_unc[sg0_bin][ThetaBin] = Esigma0;
      sg0_th[sg0_bin][ThetaBin]  = Theta;
      if(Dsigma0!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sg0_pts[sg0_bin] = ThetaBin;
    sg0_pre[sg0_bin] = Prelim;
    sg0_en[sg0_bin] = Energy;
    sg0_lo[sg0_bin] = EnergyLo;
    sg0_hi[sg0_bin] = EnergyHi;
    sg0_wt[sg0_bin] = Weight;
    sg0_sy[sg0_bin] = System;
    sg0_sc[sg0_bin] = Scale;
    strcpy(sg0_id[sg0_bin], Ident);

    //Increase energy bin counter
    sg0_bin++;
  }

  fclose(File_sg0);
  Sort_sg0(0, sg0_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sg0_bin; t++) n+=sg0_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sg0_bin; t++) if(sg0_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sg0_bin);
  for(Int_t e=0; e<sg0_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sg0_en[e], sg0_pts[e]);
    for(Int_t th=0; th<sg0_pts[e]; th++)
      printf("%f %f %f\n", sg0_th[e][th], sg0_val[e][th], sg0_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sg0(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t n0 = 0;

  for(Int_t e=0; e<sg0_bin; e++)
    if((gEnergy > sg0_lo[e]) && (gEnergy < sg0_hi[e]) && (USE_PRELIMINARY || !sg0_pre[e]))
    {
      if(bins) bins[n0] = e;
      n0++;
    }

  return n0;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sg0()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sg0 = 0.0;
  Int_t e0[EBINS];
  Int_t n0 = GetEnergyBins_sg0(e0); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigma0 data
  for(Int_t n=0; n<n0; n++) //Process all found bins
  {
    Omega = sg0_en[e0[n]];
    for(Int_t th=0; th<sg0_pts[e0[n]]; th++) //Process all data points in current bin
    {
      Theta = sg0_th[e0[n]][th];
      Meas  = sg0_sc[e0[n]]*sg0_val[e0[n]][th]*f_obs[SIG_0];
      Error = sg0_sc[e0[n]]*sg0_err[e0[n]][th]*f_obs[SIG_0];
      Theo  = sigma0(Theta, Omega);
      //printf("sg0: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sg0+=(sg0_wt[e0[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sg0;
}

//-----------------------------------------------------------------------------

void Sort_sg0(Int_t l, Int_t r) //Quicksort implementation on sg0 data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sg0_en[++i] < sg0_en[r]);
     while((sg0_en[--j] > sg0_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sg0_lo[i],  &sg0_lo[j]);
     Swap(&sg0_en[i],  &sg0_en[j]);
     Swap(&sg0_hi[i],  &sg0_hi[j]);
     Swap(&sg0_wt[i],  &sg0_wt[j]);
     Swap(&sg0_sy[i],  &sg0_sy[j]);
     Swap(&sg0_sc[i],  &sg0_sc[j]);
     Swap(&sg0_pts[i], &sg0_pts[j]);
     Swap(&sg0_pre[i], &sg0_pre[j]);
     Swap(sg0_id[i],   sg0_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sg0_val[i][n], &sg0_val[j][n]);
        Swap(&sg0_err[i][n], &sg0_err[j][n]);
        Swap(&sg0_unc[i][n], &sg0_unc[j][n]);
        Swap(&sg0_th[i][n],  &sg0_th[j][n]);
     }
   }
   Swap(&sg0_lo[i],  &sg0_lo[r]);
   Swap(&sg0_en[i],  &sg0_en[r]);
   Swap(&sg0_hi[i],  &sg0_hi[r]);
   Swap(&sg0_wt[i],  &sg0_wt[r]);
   Swap(&sg0_sy[i],  &sg0_sy[r]);
   Swap(&sg0_sc[i],  &sg0_sc[r]);
   Swap(&sg0_pts[i], &sg0_pts[r]);
   Swap(&sg0_pre[i], &sg0_pre[r]);
   Swap(sg0_id[i],   sg0_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sg0_val[i][n], &sg0_val[r][n]);
      Swap(&sg0_err[i][n], &sg0_err[r][n]);
      Swap(&sg0_unc[i][n], &sg0_unc[r][n]);
      Swap(&sg0_th[i][n],  &sg0_th[r][n]);
   }

   Sort_sg0(l, i-1);
   Sort_sg0(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sg0()
{
  Double_t Scale_sg0 = 0.0;
  Int_t e0[EBINS];
  Int_t n0 = GetEnergyBins_sg0(e0); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<n0; n++) //Process all found bins
    Scale_sg0+=(f_obs[SIG_0]-1.0)*(f_obs[SIG_0]-1.0)*sg0_pts[e0[n]]/(sg0_sy[e0[n]]*sg0_sy[e0[n]]);

  return Scale_sg0;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sg0()
{
  Int_t NPts_sg0 = 0;
  Int_t e0[EBINS];
  Int_t n0 = GetEnergyBins_sg0(e0); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<n0; n++) //Process all found bins
    NPts_sg0+=sg0_pts[e0[n]];

  return NPts_sg0;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sg0(FILE* File_sg0, Double_t* Theta, Double_t* sigma0, Double_t* Dsigma0, Double_t* Esigma0)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sg0);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigma0, Dsigma0, Esigma0);
}

//-----------------------------------------------------------------------------
