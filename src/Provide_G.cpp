#include "Provide_G.h"

//-----------------------------------------------------------------------------

void Load_G(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, G, DG, EG;
  Char_t Ident[256];
  FILE* File_G;

  printf("Loading   G  data from %s\n", Filename);
  File_G = fopen(Filename, "r");

  while(!feof(File_G))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_G, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_G, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_G(File_G, &Theta, &G, &DG, &EG)>=3)
    {
      G_val[G_bin][ThetaBin] = G;
      G_err[G_bin][ThetaBin] = DG;
      G_unc[G_bin][ThetaBin] = EG;
      G_th[G_bin][ThetaBin]  = Theta;
      if(DG!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    G_pts[G_bin] = ThetaBin;
    G_pre[G_bin] = Prelim;
    G_en[G_bin] = Energy;
    G_lo[G_bin] = EnergyLo;
    G_hi[G_bin] = EnergyHi;
    G_wt[G_bin] = Weight;
    G_sy[G_bin] = System;
    G_sc[G_bin] = Scale;
    strcpy(G_id[G_bin], Ident);

    //Increase energy bin counter
    G_bin++;
  }

  fclose(File_G);
  Sort_G(0, G_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<G_bin; t++) n+=G_pts[t];
  Int_t m = 0; for(Int_t t=0; t<G_bin; t++) if(G_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", G_bin);
  for(Int_t e=0; e<G_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, G_en[e], G_pts[e]);
    for(Int_t th=0; th<G_pts[e]; th++)
      printf("%f %f %f\n", G_th[e][th], G_val[e][th], G_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_G(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nG = 0;

  for(Int_t e=0; e<G_bin; e++)
    if((gEnergy > G_lo[e]) && (gEnergy < G_hi[e]) && (USE_PRELIMINARY || !G_pre[e]))
    {
      if(bins) bins[nG] = e;
      nG++;
    }

  return nG;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_G()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_G = 0.0;
  Int_t eG[EBINS];
  Int_t nG = GetEnergyBins_G(eG); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for G data
  for(Int_t n=0; n<nG; n++) //Process all found bins
  {
    Omega = G_en[eG[n]];
    for(Int_t th=0; th<G_pts[eG[n]]; th++) //Process all data points in current bin
    {
      Theta = G_th[eG[n]][th];
      Meas  = G_sc[eG[n]]*G_val[eG[n]][th]*f_obs[ASY_G];
      Error = G_sc[eG[n]]*G_err[eG[n]][th]*f_obs[ASY_G];
      Theo  = G(Theta, Omega);
      //printf("G: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_G+=(G_wt[eG[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_G;
}

//-----------------------------------------------------------------------------

void Sort_G(Int_t l, Int_t r) //Quicksort implementation on G data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(G_en[++i] < G_en[r]);
     while((G_en[--j] > G_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&G_lo[i],  &G_lo[j]);
     Swap(&G_en[i],  &G_en[j]);
     Swap(&G_hi[i],  &G_hi[j]);
     Swap(&G_wt[i],  &G_wt[j]);
     Swap(&G_sy[i],  &G_sy[j]);
     Swap(&G_pts[i], &G_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&G_val[i][n], &G_val[j][n]);
        Swap(&G_err[i][n], &G_err[j][n]);
        Swap(&G_unc[i][n], &G_unc[j][n]);
        Swap(&G_th[i][n],  &G_th[j][n]);
     }
   }
   Swap(&G_lo[i],  &G_lo[r]);
   Swap(&G_en[i],  &G_en[r]);
   Swap(&G_hi[i],  &G_hi[r]);
   Swap(&G_wt[i],  &G_wt[r]);
   Swap(&G_sy[i],  &G_sy[r]);
   Swap(&G_pts[i], &G_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&G_val[i][n], &G_val[r][n]);
      Swap(&G_err[i][n], &G_err[r][n]);
      Swap(&G_unc[i][n], &G_unc[r][n]);
      Swap(&G_th[i][n],  &G_th[r][n]);
   }

   Sort_G(l, i-1);
   Sort_G(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_G()
{
  Double_t Scale_G = 0.0;
  Int_t eG[EBINS];
  Int_t nG = GetEnergyBins_G(eG); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nG; n++) //Process all found bins
    Scale_G+=(1.0*G_pts[eG[n]])*(f_obs[ASY_G]-1.0)*(f_obs[ASY_G]-1.0)/(G_sy[eG[n]]*G_sy[eG[n]]);

  return Scale_G;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_G()
{
  Int_t NPts_G = 0;
  Int_t eG[EBINS];
  Int_t nG = GetEnergyBins_G(eG); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nG; n++) //Process all found bins
    NPts_G+=G_pts[eG[n]];

  return NPts_G;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_G(FILE* File_G, Double_t* Theta, Double_t* G, Double_t* DG, Double_t* EG)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_G);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, G, DG, EG);
}

//-----------------------------------------------------------------------------
