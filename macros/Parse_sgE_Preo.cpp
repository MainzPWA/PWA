static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t sgE_val[EBINS][THBINS];
Double_t sgE_err[EBINS][THBINS];
Double_t sgE_th[EBINS][THBINS];
Double_t sgE_en[EBINS];
Int_t sgE_pts[EBINS];
Int_t sgE_bin;

//-----------------------------------------------------------------------------

void Parse_sgE(Char_t* PreoFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Double_t Energy, Theta, sigmaE, DsigmaE;
  FILE* File_sgE_in;
  FILE* File_sgE_out;

  File_sgE_in = fopen(PreoFile, "r");

  sgE_bin = 0;
  while(!feof(File_sgE_in))
  {
    //Get beam energy
    if(fscanf(File_sgE_in, "E = %lf MeV\n", &Energy)!=1) break;
    ThetaBin = 0;
    while(fscanf(File_sgE_in, "%lf %lf %lf\n", &Theta, &sigmaE, &DsigmaE)==3)
    {
      sgE_val[sgE_bin][ThetaBin] = -sigmaE/2.0;
      sgE_err[sgE_bin][ThetaBin] = DsigmaE/2.0;
      sgE_th[sgE_bin][ThetaBin]  = Theta;
      if((sigmaE!=0.0) && (DsigmaE!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sgE_in);

    sgE_pts[sgE_bin] = ThetaBin;
    sgE_en[sgE_bin] = Energy;
    sgE_bin++;
  }

  fclose(File_sgE_in);

  sprintf(Buffer, "out_%s", PreoFile);
  File_sgE_out = fopen(Buffer, "w");
  for(Int_t e=0; e<sgE_bin; e++)
  {
    fprintf(File_sgE_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", sgE_en[e], Weight, System);
    for(Int_t th=0; th<sgE_pts[e]; th++)
      fprintf(File_sgE_out, "%7.3f  %10.6f  %10.6f\n", sgE_th[e][th], sgE_val[e][th], sgE_err[e][th]);
    fprintf(File_sgE_out, "-------------------------------------------\n");
  }
  fclose(File_sgE_out);
}

//-----------------------------------------------------------------------------
