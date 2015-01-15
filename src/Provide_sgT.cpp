 #include "Provide_sgT.h"

//-----------------------------------------------------------------------------

void Load_sgT(Char_t* Filename, Double_t Weight, Double_t Scale)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, sigmaT, DsigmaT, EsigmaT;
  Char_t Ident[256];
  FILE* File_sgT;

  printf("Loading sgT  data from %s\n", Filename);
  File_sgT = fopen(Filename, "r");

  while(!feof(File_sgT))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_sgT, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_sgT, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_sgT(File_sgT, &Theta, &sigmaT, &DsigmaT, &EsigmaT)>=3)
    {
      sgT_val[sgT_bin][ThetaBin] = sigmaT;
      sgT_err[sgT_bin][ThetaBin] = DsigmaT;
      sgT_unc[sgT_bin][ThetaBin] = EsigmaT;
      sgT_th[sgT_bin][ThetaBin]  = Theta;
      if(DsigmaT!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    sgT_pts[sgT_bin] = ThetaBin;
    sgT_pre[sgT_bin] = Prelim;
    sgT_en[sgT_bin] = Energy;
    sgT_lo[sgT_bin] = EnergyLo;
    sgT_hi[sgT_bin] = EnergyHi;
    sgT_wt[sgT_bin] = Weight;
    sgT_sy[sgT_bin] = System;
    sgT_sc[sgT_bin] = Scale;
    strcpy(sgT_id[sgT_bin], Ident);

    //Increase energy bin counter
    sgT_bin++;
  }

  fclose(File_sgT);
  Sort_sgT(0, sgT_bin-1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<sgT_bin; t++) n+=sgT_pts[t];
  Int_t m = 0; for(Int_t t=0; t<sgT_bin; t++) if(sgT_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", sgT_bin);
  for(Int_t e=0; e<sgT_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sgT_en[e], sgT_pts[e]);
    for(Int_t th=0; th<sgT_pts[e]; th++)
      printf("%f %f %f\n", sgT_th[e][th], sgT_val[e][th], sgT_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBins_sgT(Int_t* bins)
{
  //Build list of all energy bins covering given global energy
  Int_t nT = 0;

  for(Int_t e=0; e<sgT_bin; e++)
    if((gEnergy > sgT_lo[e]) && (gEnergy < sgT_hi[e]) && (USE_PRELIMINARY || !sgT_pre[e]))
    {
      if(bins) bins[nT] = e;
      nT++;
    }

  return nT;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_sgT()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_sgT = 0.0;
  Int_t eT[EBINS];
  Int_t nT = GetEnergyBins_sgT(eT); //Get list of all energy bins covering given global energy

  //Calculate chi^2 for sigmaT data
  for(Int_t n=0; n<nT; n++) //Process all found bins
  {
    Omega = sgT_en[eT[n]];
    for(Int_t th=0; th<sgT_pts[eT[n]]; th++) //Process all data points in current bin
    {
      Theta = sgT_th[eT[n]][th];
      Meas  = sgT_sc[eT[n]]*sgT_val[eT[n]][th]*f_obs[SIG_T];
      Error = sgT_sc[eT[n]]*sgT_err[eT[n]][th]*f_obs[SIG_T];
      Theo  = sigmaT(Theta, Omega);
      //printf("sgT: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
      ChiSq_sgT+=(sgT_wt[eT[n]]*((Meas-Theo)*(Meas-Theo)/(Error*Error)));
    }
  }
  return ChiSq_sgT;
}

//-----------------------------------------------------------------------------

void Sort_sgT(Int_t l, Int_t r) //Quicksort implementation on sgT data arrays
{
  if(r > l)
  {
    Int_t i = l-1;
    Int_t j = r;

   for(;;)
   {
     while(sgT_en[++i] < sgT_en[r]);
     while((sgT_en[--j] > sgT_en[r]) && (j>i));
     if(i>=j) break;
     Swap(&sgT_lo[i],  &sgT_lo[j]);
     Swap(&sgT_en[i],  &sgT_en[j]);
     Swap(&sgT_hi[i],  &sgT_hi[j]);
     Swap(&sgT_wt[i],  &sgT_wt[j]);
     Swap(&sgT_sy[i],  &sgT_sy[j]);
     Swap(&sgT_sc[i],  &sgT_sc[j]);
     Swap(&sgT_pts[i], &sgT_pts[j]);
     Swap(&sgT_pre[i], &sgT_pre[j]);
     Swap(sgT_id[i],   sgT_id[j]);
     for(Int_t n=0; n<THBINS; n++)
     {
        Swap(&sgT_val[i][n], &sgT_val[j][n]);
        Swap(&sgT_err[i][n], &sgT_err[j][n]);
        Swap(&sgT_unc[i][n], &sgT_unc[j][n]);
        Swap(&sgT_th[i][n],  &sgT_th[j][n]);
     }
   }
   Swap(&sgT_lo[i],  &sgT_lo[r]);
   Swap(&sgT_en[i],  &sgT_en[r]);
   Swap(&sgT_hi[i],  &sgT_hi[r]);
   Swap(&sgT_wt[i],  &sgT_wt[r]);
   Swap(&sgT_sy[i],  &sgT_sy[r]);
   Swap(&sgT_sc[i],  &sgT_sc[r]);
   Swap(&sgT_pts[i], &sgT_pts[r]);
   Swap(&sgT_pre[i], &sgT_pre[r]);
   Swap(sgT_id[i],   sgT_id[r]);
   for(Int_t n=0; n<THBINS; n++)
   {
      Swap(&sgT_val[i][n], &sgT_val[r][n]);
      Swap(&sgT_err[i][n], &sgT_err[r][n]);
      Swap(&sgT_unc[i][n], &sgT_unc[r][n]);
      Swap(&sgT_th[i][n],  &sgT_th[r][n]);
   }

   Sort_sgT(l, i-1);
   Sort_sgT(i+1, r);
  }
}

//-----------------------------------------------------------------------------

Double_t GetScale_sgT()
{
  Double_t Scale_sgT = 0.0;
  Int_t eT[EBINS];
  Int_t nT = GetEnergyBins_sgT(eT); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nT; n++) //Process all found bins
    Scale_sgT+=(f_obs[SIG_T]-1.0)*(f_obs[SIG_T]-1.0)*sgT_pts[eT[n]]/(sgT_sy[eT[n]]*sgT_sy[eT[n]]);

  return Scale_sgT;
}

//-----------------------------------------------------------------------------

Int_t GetNPts_sgT()
{
  Int_t NPts_sgT = 0;
  Int_t eT[EBINS];
  Int_t nT = GetEnergyBins_sgT(eT); //Get list of all energy bins covering given global energy

  for(Int_t n=0; n<nT; n++) //Process all found bins
    NPts_sgT+=sgT_pts[eT[n]];

  return NPts_sgT;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_sgT(FILE* File_sgT, Double_t* Theta, Double_t* sigmaT, Double_t* DsigmaT, Double_t* EsigmaT)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_sgT);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, sigmaT, DsigmaT, EsigmaT);
}

//-----------------------------------------------------------------------------
