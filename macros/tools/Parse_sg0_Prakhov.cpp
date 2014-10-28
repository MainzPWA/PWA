static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t sg0_val[EBINS][THBINS];
Double_t sg0_err[EBINS][THBINS];
Double_t sg0_unc[EBINS][THBINS];
Double_t sg0_th[EBINS][THBINS];
Double_t sg0_lo[EBINS];
Double_t sg0_hi[EBINS];
Double_t sg0_en[EBINS];
Double_t sg0_sy[EBINS];
Int_t sg0_pts[EBINS];
Int_t sg0_bin;

//-----------------------------------------------------------------------------

void Parse_sg0_v2(Char_t* PrakhovFile, Char_t* ID)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Double_t Energy, Theta, sigma0, Dsigma0, Esigma0;
  Double_t Dummy1, Dummy2, Dummy3, Dummy4;
  Double_t EnergyLo, EnergyHi;
  FILE* File_sg0_in;
  FILE* File_sg0_out;

  File_sg0_in = fopen(PrakhovFile, "r");

  sg0_bin = 0;
  while(!feof(File_sg0_in))
  {
    //Get lo/hi beam energy
    if(fscanf(File_sg0_in, "Eg=%lf-%lf W=%lf-%lf\n", &EnergyLo, &EnergyHi, &Dummy3, &Dummy4)!=4) break;
    //Get mean beam energy
    if(fscanf(File_sg0_in, "Eg=%lf+/-%lf W=%lf+/-%lf\n", &Energy, &Dummy1, &Dummy2, &Dummy3)!=4) break;
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sg0_in);
    ThetaBin = 0;
    while(fscanf(File_sg0_in, "%lf %lf %lf %lf %lf\n", &Dummy1, &Theta, &sigma0, &Dsigma0, &Esigma0)==5)
    {
      sg0_val[sg0_bin][ThetaBin] = sigma0;
      sg0_err[sg0_bin][ThetaBin] = Dsigma0;
      sg0_unc[sg0_bin][ThetaBin] = TMath::Sqrt(Esigma0*Esigma0 - Dsigma0*Dsigma0);
      sg0_th[sg0_bin][ThetaBin]  = Theta;
      if(Dsigma0!=0.0) ThetaBin++; //Accept only 'existing' data points
    }

    for(Int_t t=0; t<ThetaBin; t++)
      sg0_sy[sg0_bin]+=sg0_unc[sg0_bin][t];
    sg0_sy[sg0_bin]/=ThetaBin;

    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sg0_in);
    fgets(Buffer, sizeof(Buffer), File_sg0_in);

    sg0_pts[sg0_bin] = ThetaBin;
    sg0_en[sg0_bin] = Energy;
    sg0_lo[sg0_bin] = EnergyLo;
    sg0_hi[sg0_bin] = EnergyHi;
    sg0_bin++;
  }
  fclose(File_sg0_in);

  sprintf(Buffer, "out_%s", PrakhovFile);
  File_sg0_out = fopen(Buffer, "w");
  for(Int_t e=0; e<sg0_bin; e++)
  {
    fprintf(File_sg0_out, "E = %8.3f MeV, E_lo = %8.3f MeV, E_hi = %8.3f MeV\n", sg0_en[e], sg0_lo[e], sg0_hi[e]);
    fprintf(File_sg0_out, "Systematic = %6.4f, Preliminary = 0, %s\n", sg0_sy[e], ID);

    for(Int_t th=0; th<sg0_pts[e]; th++)
      fprintf(File_sg0_out, "%7.3f   %10.6f   %10.6f   %10.6f\n", sg0_th[e][th], sg0_val[e][th], sg0_err[e][th], sg0_unc[e][th]);
    fprintf(File_sg0_out, "----------------------------------------------------------\n");
  }
  fclose(File_sg0_out);
}

//-----------------------------------------------------------------------------

void Parse_sg0_v1(Char_t* PrakhovFile, Char_t* ID)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Double_t Energy, Theta, sigma0, Dsigma0;
  Double_t Dummy1, Dummy2, Dummy3, Dummy4;
  Double_t EnergyLo, EnergyHi;
  FILE* File_sg0_in;
  FILE* File_sg0_out;

  File_sg0_in = fopen(PrakhovFile, "r");

  sg0_bin = 0;
  while(!feof(File_sg0_in))
  {
    //Get lo/hi beam energy
    if(fscanf(File_sg0_in, "Eg=%lf-%lf W=%lf-%lf\n", &EnergyLo, &EnergyHi, &Dummy3, &Dummy4)!=4) break;
    //Get mean beam energy
    if(fscanf(File_sg0_in, "Eg=%lf+/-%lf W=%lf+/-%lf\n", &Energy, &Dummy1, &Dummy2, &Dummy3)!=4) break;
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sg0_in);
    ThetaBin = 0;
    while(fscanf(File_sg0_in, "%lf %lf %lf %lf\n", &Dummy1, &Theta, &sigma0, &Dsigma0)==4)
    {
      sg0_val[sg0_bin][ThetaBin] = sigma0;
      sg0_err[sg0_bin][ThetaBin] = Dsigma0;
      sg0_th[sg0_bin][ThetaBin]  = Theta;
      if(Dsigma0!=0.0) ThetaBin++; //Accept only 'existing' data points
    }

    sg0_sy[sg0_bin] = 0.05;

    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sg0_in);
    fgets(Buffer, sizeof(Buffer), File_sg0_in);

    sg0_pts[sg0_bin] = ThetaBin;
    sg0_en[sg0_bin] = Energy;
    sg0_lo[sg0_bin] = EnergyLo;
    sg0_hi[sg0_bin] = EnergyHi;
    sg0_bin++;
  }
  fclose(File_sg0_in);

  sprintf(Buffer, "out_%s", PrakhovFile);
  File_sg0_out = fopen(Buffer, "w");
  for(Int_t e=0; e<sg0_bin; e++)
  {
    fprintf(File_sg0_out, "E = %8.3f MeV, E_lo = %8.3f MeV, E_hi = %8.3f MeV\n", sg0_en[e], sg0_lo[e], sg0_hi[e]);
    fprintf(File_sg0_out, "Systematic = %6.4f, Preliminary = 0, %s\n", sg0_sy[e], ID);

    for(Int_t th=0; th<sg0_pts[e]; th++)
      fprintf(File_sg0_out, "%7.3f   %10.6f   %10.6f\n", sg0_th[e][th], sg0_val[e][th], sg0_err[e][th]);
    fprintf(File_sg0_out, "----------------------------------------------------------\n");
  }
  fclose(File_sg0_out);
}

//-----------------------------------------------------------------------------
