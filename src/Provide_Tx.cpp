#include "Provide_Tx.h"

//-----------------------------------------------------------------------------

void Load_Tx(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Tx, DTx, ETx;
  Char_t Ident[256];
  FILE* File_Tx;

  printf("Loading   Tx data from %s\n", Filename);
  File_Tx = fopen(Filename, "r");

  while(!feof(File_Tx))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Tx, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Tx, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Tx(File_Tx, &Theta, &Tx, &DTx, &ETx)>=3)
    {
      Tx_val[Tx_bin][ThetaBin] = Tx;
      Tx_err[Tx_bin][ThetaBin] = DTx;
      Tx_unc[Tx_bin][ThetaBin] = ETx;
      Tx_th[Tx_bin][ThetaBin]  = Theta;
      if(DTx!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Tx_pts[Tx_bin] = ThetaBin;
    Tx_pre[Tx_bin] = Prelim;
    Tx_en[Tx_bin] = Energy;
    Tx_lo[Tx_bin] = EnergyLo;
    Tx_hi[Tx_bin] = EnergyHi;
    Tx_wt[Tx_bin] = Weight;
    Tx_sy[Tx_bin] = System;
    Tx_sc[Tx_bin] = Scale;
    strcpy(Tx_id[Tx_bin], Ident);

    //Increase energy bin counter
    Tx_bin++;
  }

  fclose(File_Tx);
  Sort_Tx(0, Tx_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Tx_bin; t++) n+=Tx_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Tx_bin; t++) if(Tx_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Tx_bin);
  for(Int_t e=0; e<Tx_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Tx_en[e], Tx_pts[e]);
    for(Int_t th=0; th<Tx_pts[e]; th++)
      printf("%f %f %f\n", Tx_th[e][th], Tx_val[e][th], Tx_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_Tx(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nTx = 0;

  for(Int_t e=0; e<Tx_bin; e++)
    if((gEnergy > Tx_lo[e]) && (gEnergy < Tx_hi[e]) && (USE_PRELIMINARY || !Tx_pre[e]))
    {
      if(bins) bins[nTx] = e;
      nTx++;
    }

  return nTx;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_Tx()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_Tx = 0.0;
  Int_t eTx[EBINS];
  Int_t nTx = GetEnergyBins_Tx(eTx); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for Tx data
  for(Int_t n=0; n<nTx; n++) //Process all found bins
  {
    Omega = Tx_en[eTx[n]];
    for(Int_t th=0; th<Tx_pts[eTx[n]]; th++) //Process all data points in current bin
    {
      Theta = Tx_th[eTx[n]][th];
      Meas  = Tx_sc[eTx[n]]*Tx_val[eTx[n]][th]*f_obs[ASY_TX];
      Error = Tx_sc[eTx[n]]*Tx_err[eTx[n]][th]*f_obs[ASY_TX];
      Theo  = Tx(Theta, Omega);
      //printf("Tx: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_Tx+=(Tx_wt[eTx[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_Tx;
}

//-----------------------------------------------------------------------------

void Sort_Tx(Int_t l, Int_t r) //Quicksort implementation on Tx data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(Tx_en[++i] < Tx_en[r]);
     while((Tx_en[--j] > Tx_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&Tx_lo[i],  &Tx_lo[j]);
     Swap(&Tx_en[i],  &Tx_en[j]);
     Swap(&Tx_hi[i],  &Tx_hi[j]);
     Swap(&Tx_wt[i],  &Tx_wt[j]);
     Swap(&Tx_sy[i],  &Tx_sy[j]);
     Swap(&Tx_pts[i], &Tx_pts[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&Tx_val[i][n], &Tx_val[j][n]);
        Swap(&Tx_err[i][n], &Tx_err[j][n]);
        Swap(&Tx_unc[i][n], &Tx_unc[j][n]);
        Swap(&Tx_th[i][n],  &Tx_th[j][n]);
     }
   }
   Swap(&Tx_lo[i],  &Tx_lo[r]);
   Swap(&Tx_en[i],  &Tx_en[r]);
   Swap(&Tx_hi[i],  &Tx_hi[r]);
   Swap(&Tx_wt[i],  &Tx_wt[r]);
   Swap(&Tx_sy[i],  &Tx_sy[r]);
   Swap(&Tx_pts[i], &Tx_pts[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&Tx_val[i][n], &Tx_val[r][n]);
      Swap(&Tx_err[i][n], &Tx_err[r][n]);
      Swap(&Tx_unc[i][n], &Tx_unc[r][n]);
      Swap(&Tx_th[i][n],  &Tx_th[r][n]);
   }

   Sort_Tx(l, i-1);
   Sort_Tx(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_Tx()
{
  Double_t Scale_Tx = 0.0;
  Int_t eTx[EBINS];
  Int_t nTx = GetEnergyBins_Tx(eTx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nTx; n++) //Process all found bins
    Scale_Tx+=(f_obs[ASY_TX]-1.0)*(f_obs[ASY_TX]-1.0)*Tx_pts[eTx[n]]/(Tx_sy[eTx[n]]*Tx_sy[eTx[n]]);

  return Scale_Tx;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_Tx()
{
  Int_t NPts_Tx = 0;
  Int_t eTx[EBINS];
  Int_t nTx = GetEnergyBins_Tx(eTx); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nTx; n++) //Process all found bins
    NPts_Tx+=Tx_pts[eTx[n]];

  return NPts_Tx;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Tx(FILE* File_Tx, Double_t* Theta, Double_t* Tx, Double_t* DTx, Double_t* ETx)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Tx);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Tx, DTx, ETx);
}

//-----------------------------------------------------------------------------
