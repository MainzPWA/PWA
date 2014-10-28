static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t F_val[EBINS][THBINS];
Double_t F_err[EBINS][THBINS];
Double_t F_th[EBINS][THBINS];
Double_t F_en[EBINS];
Int_t F_pts[EBINS];
Int_t F_bin;

//-----------------------------------------------------------------------------

void Parse_F(Char_t* KashevFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Double_t Energy, CosTheta, F, DF;
  FILE* File_F_in;
  FILE* File_F_out;

  File_F_in = fopen(KashevFile, "r");

  F_bin = 0;
  while(!feof(File_F_in))
  {
    //Get beam energy
    if(fscanf(File_F_in, "#       E_gamma = %lf MeV\n", &Energy)!=1) break;
    ThetaBin = 0;
    while(fscanf(File_F_in, "%lf %lf %lf\n", &CosTheta, &F, &DF)==3)
    {
      F_val[F_bin][ThetaBin] = F;
      F_err[F_bin][ThetaBin] = DF;
      F_th[F_bin][ThetaBin]  = TMath::ACos(CosTheta)*TMath::RadToDeg();
      if(DF!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_F_in);

    F_pts[F_bin] = ThetaBin;
    F_en[F_bin] = Energy;
    F_bin++;
  }

  fclose(File_F_in);

  sprintf(Buffer, "out_%s", KashevFile);
  File_F_out = fopen(Buffer, "w");
  for(Int_t e=0; e<F_bin; e++)
  {
    fprintf(File_F_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", F_en[e], Weight, System);
    for(Int_t th=0; th<F_pts[e]; th++)
      fprintf(File_F_out, "%7.3f  %10.6f  %10.6f\n", F_th[e][th], F_val[e][th], F_err[e][th]);
      fprintf(File_F_out, "-------------------------------------------\n");
  }
  fclose(File_F_out);
}

//-----------------------------------------------------------------------------
