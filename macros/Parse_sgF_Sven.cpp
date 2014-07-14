static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t sgF_val[EBINS][THBINS];
Double_t sgF_err[EBINS][THBINS];
Double_t sgF_th[EBINS][THBINS];
Double_t sgF_en[EBINS];
Int_t sgF_pts[EBINS];
Int_t sgF_bin;

//-----------------------------------------------------------------------------

void Parse_sgF(Char_t* SvenFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Int_t Dummy1, Dummy2;
  Double_t Energy, Theta, DTheta, sigmaF, DsigmaF;
  Double_t Dummy3;
  FILE* File_sgF_in;
  FILE* File_sgF_out;

  File_sgF_in = fopen(SvenFile, "r");

  sgF_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), File_sgF_in);
  fgets(Buffer, sizeof(Buffer), File_sgF_in);
  while(!feof(File_sgF_in))
  {
    //Skip 3 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sgF_in);
    fgets(Buffer, sizeof(Buffer), File_sgF_in);
    fgets(Buffer, sizeof(Buffer), File_sgF_in);
    //Get beam energy
    if(fscanf(File_sgF_in, "Channel %d to %d: %lf +- %lf MeV (preliminary)\n", &Dummy1, &Dummy2, &Energy, &Dummy3)!=4) break;
    //Skip 1 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sgF_in);
    ThetaBin = 0;
    while(fscanf(File_sgF_in, "%lf %lf %lf %lf\n", &Theta, &DTheta, &sigmaF, &DsigmaF)==4)
    {
      sgF_val[sgF_bin][ThetaBin] = sigmaF;
      sgF_err[sgF_bin][ThetaBin] = DsigmaF;
      sgF_th[sgF_bin][ThetaBin]  = Theta;
      if((sigmaF!=0.0) && (DsigmaF!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 5 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sgF_in);
    fgets(Buffer, sizeof(Buffer), File_sgF_in);
    fgets(Buffer, sizeof(Buffer), File_sgF_in);
    fgets(Buffer, sizeof(Buffer), File_sgF_in);

    sgF_pts[sgF_bin] = ThetaBin;
    sgF_en[sgF_bin] = Energy;
    sgF_bin++;
  }

  fclose(File_sgF_in);

  sprintf(Buffer, "out_%s", SvenFile);
  File_sgF_out = fopen(Buffer, "w");
  for(Int_t e=0; e<sgF_bin; e++)
  {
    fprintf(File_sgF_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", sgF_en[e], Weight, System);
    for(Int_t th=0; th<sgF_pts[e]; th++)
      fprintf(File_sgF_out, "%7.3f  %10.6f  %10.6f\n", sgF_th[e][th], sgF_val[e][th], sgF_err[e][th]);
    fprintf(File_sgF_out, "-------------------------------------------\n");
  }
  fclose(File_sgF_out);
}

//-----------------------------------------------------------------------------
