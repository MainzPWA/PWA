
static const Int_t EBINS=1024;
static const Int_t THBINS=256;

Double_t Obs1_val[EBINS][THBINS];
Double_t Obs1_err[EBINS][THBINS];
Double_t Obs1_unc[EBINS][THBINS];
Double_t Obs1_th[EBINS][THBINS];
Double_t Obs1_lo[EBINS];
Double_t Obs1_en[EBINS];
Double_t Obs1_hi[EBINS];
Double_t Obs1_sy[EBINS];
Int_t Obs1_pre[EBINS];
Int_t Obs1_pts[EBINS];
Int_t Obs1_bin;
Char_t Obs1_id[EBINS][256];

Double_t Obs2_val[EBINS][THBINS];
Double_t Obs2_err[EBINS][THBINS];
Double_t Obs2_unc[EBINS][THBINS];
Double_t Obs2_th[EBINS][THBINS];
Double_t Obs2_lo[EBINS];
Double_t Obs2_en[EBINS];
Double_t Obs2_hi[EBINS];
Double_t Obs2_sy[EBINS];
Int_t Obs2_pre[EBINS];
Int_t Obs2_pts[EBINS];
Int_t Obs2_bin;
Char_t Obs2_id[EBINS][256];

//-----------------------------------------------------------------------------

void Load1(Char_t* Filename)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Obs, DObs, EObs;
  Char_t Ident[256];
  FILE* File_Obs1;

  printf("Loading observable data from %s\n", Filename);
  File_Obs1 = fopen(Filename, "r");

  Obs1_bin = 0;
  while(!feof(File_Obs1))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Obs1, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Obs1, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Obs(File_Obs1, &Theta, &Obs, &DObs, &EObs)>=3)
    {
      Obs1_val[Obs1_bin][ThetaBin] = Obs;
      Obs1_err[Obs1_bin][ThetaBin] = DObs;
      Obs1_unc[Obs1_bin][ThetaBin] = EObs;
      Obs1_th[Obs1_bin][ThetaBin]  = Theta;
      if(DObs!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Obs1_pts[Obs1_bin] = ThetaBin;
    Obs1_pre[Obs1_bin] = Prelim;
    Obs1_en[Obs1_bin] = Energy;
    Obs1_lo[Obs1_bin] = EnergyLo;
    Obs1_hi[Obs1_bin] = EnergyHi;
    Obs1_sy[Obs1_bin] = System;
    strcpy(Obs1_id[Obs1_bin], Ident);

    //Increase energy bin counter
    Obs1_bin++;
  }

  fclose(File_Obs1);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Obs1_bin; t++) n+=Obs1_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Obs1_bin; t++) if(Obs1_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Obs1_bin);
  for(Int_t e=0; e<Obs1_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Obs1_en[e], Obs1_pts[e]);
    for(Int_t th=0; th<Obs1_pts[e]; th++)
      printf("%f %f %f\n", Obs1_th[e][th], Obs1_val[e][th], Obs1_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

void Load2(Char_t* Filename)
{
  Int_t Prelim;
  Int_t ThetaBin;
  Double_t Energy, EnergyLo, EnergyHi;
  Double_t System;
  Double_t Theta, Obs, DObs, EObs;
  Char_t Ident[256];
  FILE* File_Obs2;

  printf("Loading observable data from %s\n", Filename);
  File_Obs2 = fopen(Filename, "r");

  Obs2_bin = 0;
  while(!feof(File_Obs2))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Obs2, "E = %lf MeV, E_lo = %lf MeV, E_hi = %lf MeV\n", &Energy, &EnergyLo, &EnergyHi)!=3) break;
    if(fscanf(File_Obs2, "Systematic = %lf, Preliminary = %d, %s\n", &System, &Prelim, Ident)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Obs(File_Obs2, &Theta, &Obs, &DObs, &EObs)>=3)
    {
      Obs2_val[Obs2_bin][ThetaBin] = Obs;
      Obs2_err[Obs2_bin][ThetaBin] = DObs;
      Obs2_unc[Obs2_bin][ThetaBin] = EObs;
      Obs2_th[Obs2_bin][ThetaBin]  = Theta;
      if(DObs!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    //Store data for this energy bin
    Obs2_pts[Obs2_bin] = ThetaBin;
    Obs2_pre[Obs2_bin] = Prelim;
    Obs2_en[Obs2_bin] = Energy;
    Obs2_lo[Obs2_bin] = EnergyLo;
    Obs2_hi[Obs2_bin] = EnergyHi;
    Obs2_sy[Obs2_bin] = System;
    strcpy(Obs2_id[Obs2_bin], Ident);

    //Increase energy bin counter
    Obs2_bin++;
  }

  fclose(File_Obs2);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Obs2_bin; t++) n+=Obs2_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Obs2_bin; t++) if(Obs2_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);
  return;

  //Debug output
  printf("EBins: %d\n", Obs2_bin);
  for(Int_t e=0; e<Obs2_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Obs2_en[e], Obs2_pts[e]);
    for(Int_t th=0; th<Obs2_pts[e]; th++)
      printf("%f %f %f\n", Obs2_th[e][th], Obs2_val[e][th], Obs2_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

Int_t ReadLine_Obs(FILE* File_Obs, Double_t* Theta, Double_t* Obs, Double_t* DObs, Double_t* EObs)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Obs);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Obs, DObs, EObs);
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin1(Double_t Energy)
{
  //Get energy bin for sigma0 for given global energy
  Double_t Min = 1e38;
  Int_t e1 = 0;

  for(Int_t e=0; e<Obs1_bin; e++)
    if(fabs(Obs1_en[e] - Energy) < Min)
    {
      Min = fabs(Obs1_en[e] - Energy);
      e1 = e;
    }
  return e1;
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin2(Double_t Energy)
{
  //Get energy bin for sigma0 for given global energy
  Double_t Min = 1e38;
  Int_t e2 = 0;

  for(Int_t e=0; e<Obs2_bin; e++)
    if(fabs(Obs2_en[e] - Energy) < Min)
    {
      Min = fabs(Obs2_en[e] - Energy);
      e2 = e;
    }
  return e2;
}

//-----------------------------------------------------------------------------

void Compare(Char_t* Obs, Double_t Energy)
{
  Char_t Buffer[256];
  TH1F* Hist_Sven   = new TH1F("Sven",   "Sven",   18, 0.0, 180.0);
  TH1F* Hist_Jannes = new TH1F("Jannes", "Jannes", 18, 0.0, 180.0);

  sprintf(Buffer, "Sven_fine/%s.txt", Obs);
  Load1(Buffer);
  sprintf(Buffer, "data/%s.txt", Obs);
  Load2(Buffer);

  Int_t Bin_Sven   = GetEnergyBin1(Energy);
  Int_t Bin_Jannes = GetEnergyBin2(Energy);

  for(Int_t t=0; t<18; t++)
  {
    Hist_Sven->SetBinContent(t+1, Obs1_val[Bin_Sven][t]);
    Hist_Sven->SetBinError(t+1, Obs1_err[Bin_Sven][t]);
    Hist_Jannes->SetBinContent(t+1, Obs2_val[Bin_Jannes][t]);
    Hist_Jannes->SetBinError(t+1, Obs2_err[Bin_Jannes][t]);
  }
  Hist_Sven->SetLineColor(kBlack);
  Hist_Sven->SetMarkerColor(kBlack);
  Hist_Sven->SetMarkerStyle(20);
  Hist_Sven->SetMarkerSize(0.8);
  Hist_Sven->Draw();

  Hist_Jannes->SetLineColor(kRed);
  Hist_Jannes->SetMarkerColor(kRed);
  Hist_Jannes->SetMarkerStyle(20);
  Hist_Jannes->SetMarkerSize(0.8);
  Hist_Jannes->Draw("same");
}

//-----------------------------------------------------------------------------
