#include "Parse_G.h"

//-----------------------------------------------------------------------------

void Parse_G()
{
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, G, DG;
  FILE* File_G;

  printf("Loading   G  data... ");
  File_G = fopen("data/G.txt", "r");

  G_bin = 0;
  while(!feof(File_G))
  {
    //Get beam energy
    if(fscanf(File_G, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    //Check if this energy already exists
    EnergyBin = ExistEnergyBin_G(Energy);
    if(EnergyBin==G_bin) //This energy is new (index is at end)...
      ThetaBin = 0; //...hence we will set up a new energy bin
    else //This energy already exists...
      ThetaBin = G_pts[EnergyBin]; //..hence we append to the existing energy bin

    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_G(File_G, &Theta, &G, &DG)==3)
    {
      G_val[EnergyBin][ThetaBin] = G;
      G_err[EnergyBin][ThetaBin] = DG;
      G_th[EnergyBin][ThetaBin]  = Theta;
      if(DG!=0.0) ThetaBin++; //Accept only 'existing' data points
    }

    G_pts[EnergyBin] = ThetaBin;
    G_en[EnergyBin] = Energy;
    G_wt[EnergyBin] = Weight;
    G_sy[EnergyBin] = System;

    //Increase energy bin counter, if energy bin is newly set up
    if(EnergyBin==G_bin)
      G_bin++;
  }

  fclose(File_G);
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

Int_t GetEnergyBin_G()
{
  //Get energy bin for G for given global energy
  Double_t Min = 1e38;
  Int_t eG = 0;

  for(Int_t e=0; e<G_bin; e++)
    if(fabs(G_en[e] - gEnergy) < Min)
    {
      Min = fabs(G_en[e] - gEnergy);
      eG = e;
    }
  return eG;
}

//-----------------------------------------------------------------------------

Int_t ExistEnergyBin_G(Double_t Energy)
{
  //Check for existing energy bin
  Int_t eG = G_bin;

  for(Int_t e=0; e<G_bin; e++)
    if(G_en[e]==Energy)
      eG = e;

  return eG;
}

//-----------------------------------------------------------------------------

Double_t GetChiSq_G()
{
  Double_t Meas;
  Double_t Theo;
  Double_t Error;
  Double_t Theta;
  Double_t Omega;
  Double_t ChiSq_G;
  Int_t eG = GetEnergyBin_G();

  //Calculate chi^2 for G data
  ChiSq_G = 0.0;
  Omega = G_en[eG];
  for(Int_t th=0; th<G_pts[eG]; th++)
  {
    Theta = G_th[eG][th];
    Meas  = G_val[eG][th]*f_obs[ASY_G];
    Error = G_err[eG][th]*f_obs[ASY_G];
    Theo  = G(Theta, Omega);
    //printf("G: %f: %f %f  = %f\n", Theta, Theo, Meas, Theo/Meas);
    ChiSq_G+=((Meas-Theo)*(Meas-Theo)/(Error*Error));
  }

  return G_wt[eG]*ChiSq_G; //Apply weight factor for each dataset
}

//-----------------------------------------------------------------------------

Double_t GetScale_G()
{
  Int_t eG = GetEnergyBin_G();
  Double_t Scale_G = (1.0*G_pts[eG])*(f_obs[ASY_G]-1.0)*(f_obs[ASY_G]-1.0)/(G_sy[eG]*G_sy[eG]);
  return Scale_G;
}

//-----------------------------------------------------------------------------

Int_t ReadLine_G(FILE* File_G, Double_t* Theta, Double_t* G, Double_t* DG)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_G);
  return sscanf(Buffer, "%lf %lf %lf", Theta, G, DG);
}

//-----------------------------------------------------------------------------
