
static const Int_t EBINS=512;
static const Int_t THBINS=256;

Double_t Obs_val[EBINS][THBINS];
Double_t Obs_err[EBINS][THBINS];
Double_t Obs_unc[EBINS][THBINS];
Double_t Obs_th[EBINS][THBINS];
Double_t Obs_lo[EBINS];
Double_t Obs_en[EBINS];
Double_t Obs_hi[EBINS];
Double_t Obs_wt[EBINS];
Double_t Obs_sy[EBINS];
Double_t Obs_sc[EBINS];
Int_t Obs_pre[EBINS];
Int_t Obs_pts[EBINS];
Int_t Obs_bin;
Char_t Obs_id[EBINS][256];


void Convert(Char_t* Filename, Char_t* ID)
{
  Char_t Buff[256];
  Double_t E[512];
  Double_t DE[512];
  Double_t Val, Err;
  Int_t Nch = 0;
  Double_t En, En_lo, En_hi;
  Double_t Wt, Sy;
  Double_t th, X, DX, EX;
  Int_t ThetaBin;

  //FILE* File_Eg = fopen("ChanEn450.txt", "r");
  //FILE* File_Eg = fopen("ChanEn883.txt", "r");
  //FILE* File_Eg = fopen("ChanEn1508.txt", "r");
  FILE* File_Eg = fopen("ChanEn.txt", "r");
  while(!feof(File_Eg))
  {
    if(fscanf(File_Eg, "%lf %lf\n", &Val, &Err)==2)
    {
      E[Nch]  = Val;
      DE[Nch] = Err;
      Nch++;
    }
  }

  Obs_bin = 0;
  FILE* File_Obs = fopen(Filename, "r");

  while(!feof(File_Obs))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Obs, "E =  %lf MeV, Wght = %lf, Syst = %lf\n", &En, &Wt, &Sy)!=3) break;

    ThetaBin = 0;
    //This will read lines from file until end-of-entry marker (e.g. "---...---" line) is found
    while(ReadLine_Obs(File_Obs, &th, &X, &DX, &EX)>=3)
    {
      Obs_val[Obs_bin][ThetaBin] = X;
      Obs_err[Obs_bin][ThetaBin] = DX;
      Obs_unc[Obs_bin][ThetaBin] = EX;
      Obs_th[Obs_bin][ThetaBin]  = th;
      if(DX!=0.0) ThetaBin++; //Accept only 'existing' data points (i.e. with finite error)
    }

    for(Int_t e=0; e<Nch; e++)
      if(fabs(En-E[e]) < 0.01)
      {
        En_lo = En - (DE[e]/2.0);
        En_hi = En + (DE[e]/2.0);
        //En_lo = En - (DE[e]/1.0);
        //En_hi = En + (DE[e]/1.0);
      }

    //Store data for this energy bin
    Obs_pts[Obs_bin] = ThetaBin;
    Obs_en[Obs_bin] = En;
    Obs_lo[Obs_bin] = En_lo;
    Obs_hi[Obs_bin] = En_hi;
    Obs_sy[Obs_bin] = Sy;

    //Increase energy bin counter
    Obs_bin++;
  }

  fclose(File_Obs);

  sprintf(Buff, "new_%s", Filename);
  FILE* File_New = fopen(Buff, "w");

  for(Int_t e=0; e<Obs_bin; e++)
  {
    fprintf(File_New, "E = %8.3f MeV, E_lo = %8.3f MeV, E_hi = %8.3f MeV\n", Obs_en[e], Obs_lo[e], Obs_hi[e]);
    fprintf(File_New, "Systematic = %8.6f, Preliminary = 0, %s\n", Obs_sy[e], ID);
    for(Int_t t=0; t<Obs_pts[e]; t++)
      //fprintf(File_New, "%7.3f   %10.6f   %10.6f   %10.6f\n", Obs_th[e][t], Obs_val[e][t], Obs_err[e][t], Obs_unc[e][t]);
      fprintf(File_New, "%7.3f   %10.6f   %10.6f\n", Obs_th[e][t], Obs_val[e][t], Obs_err[e][t]);
    sprintf(Buff, "Systematic = %8.6f, %s\n", Obs_sy[e], ID);
    fprintf(File_New, "----------------------------------------------------------\n");
  }
  fclose(File_New);

}

//-----------------------------------------------------------------------------

Int_t ReadLine_Obs(FILE* File_Obs, Double_t* Theta, Double_t* Obs, Double_t* DObs, Double_t* EObs)
{
  Char_t Buffer[1024];

  fgets(Buffer, sizeof(Buffer), File_Obs);
  return sscanf(Buffer, "%lf %lf %lf %lf", Theta, Obs, DObs, EObs);
}

//-----------------------------------------------------------------------------

void Convert_Tuzla(Char_t* Filename, Char_t* ID)
{
  Char_t Buff[256];
  Double_t E[512];
  Double_t DE[512];
  Double_t Val, Err;
  Int_t Nch = 0;
  Double_t En, En_lo, En_hi;
  Double_t Wt, Sy;
  Double_t th, X, DX, EX;
  Int_t ThetaBin;

  //FILE* File_Eg = fopen("ChanEn450.txt", "r");
  //FILE* File_Eg = fopen("ChanEn883.txt", "r");
  //FILE* File_Eg = fopen("ChanEn1508.txt", "r");
  FILE* File_Eg = fopen("ChanEn.txt", "r");
  while(!feof(File_Eg))
  {
    if(fscanf(File_Eg, "%lf %lf\n", &Val, &Err)==2)
    {
      E[Nch]  = Val;
      DE[Nch] = Err;
      Nch++;
    }
  }

  Obs_bin = 0;
  FILE* File_Obs = fopen(Filename, "r");

  while(!feof(File_Obs))
  {
    //Get header informations (energy, weight, ID, ...)
    if(fscanf(File_Obs, "%lf %d\n", &En, &ThetaBin)!=2) break;
printf("%f %d\n", En, ThetaBin);

    for(Int_t n=0; n<ThetaBin; n++)
    {
      if(fscanf(File_Obs, "%lf %lf %lf\n", &th, &X, &DX)!=3) break;
printf("%f %f %f\n", th, X, DX);

      Obs_val[Obs_bin][n] = X;
      Obs_err[Obs_bin][n] = DX;
      Obs_unc[Obs_bin][n] = 0.0;
      Obs_th[Obs_bin][n]  = TMath::ACos(th)*TMath::RadToDeg();
    }

    for(Int_t e=0; e<Nch; e++)
      if(fabs(En-E[e]) < 0.01)
      {
        En_lo = En - (DE[e]/2.0);
        En_hi = En + (DE[e]/2.0);
        //En_lo = En - (DE[e]/1.0);
        //En_hi = En + (DE[e]/1.0);
      }

    //Store data for this energy bin
    Obs_pts[Obs_bin] = ThetaBin;
    Obs_en[Obs_bin] = En;
    Obs_lo[Obs_bin] = En_lo;
    Obs_hi[Obs_bin] = En_hi;
    Obs_sy[Obs_bin] = 0.05;

    //Increase energy bin counter
    Obs_bin++;
  }

  fclose(File_Obs);
printf("%d\n", Obs_bin);

  sprintf(Buff, "new_%s", Filename);
  FILE* File_New = fopen(Buff, "w");

  for(Int_t e=0; e<Obs_bin; e++)
  {
    fprintf(File_New, "E = %8.3f MeV, E_lo = %8.3f MeV, E_hi = %8.3f MeV\n", Obs_en[e], Obs_lo[e], Obs_hi[e]);
    fprintf(File_New, "Systematic = %8.6f, Preliminary = 0, %s\n", Obs_sy[e], ID);
    for(Int_t t=0; t<Obs_pts[e]; t++)
      //fprintf(File_New, "%7.3f   %10.6f   %10.6f   %10.6f\n", Obs_th[e][t], Obs_val[e][t], Obs_err[e][t], Obs_unc[e][t]);
      fprintf(File_New, "%7.3f   %10.6f   %10.6f\n", Obs_th[e][t], Obs_val[e][t], Obs_err[e][t]);
    sprintf(Buff, "Systematic = %8.6f, %s\n", Obs_sy[e], ID);
    fprintf(File_New, "----------------------------------------------------------\n");
  }
  fclose(File_New);

}

//-----------------------------------------------------------------------------
