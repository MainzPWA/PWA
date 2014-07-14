static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t sgT_val[EBINS][THBINS];
Double_t sgT_err[EBINS][THBINS];
Double_t sgT_th[EBINS][THBINS];
Double_t sgT_en[EBINS];
Int_t sgT_pts[EBINS];
Int_t sgT_bin;

//-----------------------------------------------------------------------------

void Parse_sgT(Char_t* SvenFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Int_t Dummy1, Dummy2;
  Double_t Energy, Theta, DTheta, sigmaT, DsigmaT;
  Double_t Dummy3;
  FILE* File_sgT_in;
  FILE* File_sgT_out;

  File_sgT_in = fopen(SvenFile, "r");

  sgT_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), File_sgT_in);
  fgets(Buffer, sizeof(Buffer), File_sgT_in);
  while(!feof(File_sgT_in))
  {
    //Skip 3 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sgT_in);
    fgets(Buffer, sizeof(Buffer), File_sgT_in);
    fgets(Buffer, sizeof(Buffer), File_sgT_in);
    //Get beam energy
    if(fscanf(File_sgT_in, "Channel %d to %d: %lf +- %lf MeV (preliminary)\n", &Dummy1, &Dummy2, &Energy, &Dummy3)!=4) break;
    //Skip 1 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sgT_in);
    ThetaBin = 0;
    while(fscanf(File_sgT_in, "%lf %lf %lf %lf\n", &Theta, &DTheta, &sigmaT, &DsigmaT)==4)
    {
      sgT_val[sgT_bin][ThetaBin] = sigmaT;
      sgT_err[sgT_bin][ThetaBin] = DsigmaT;
      sgT_th[sgT_bin][ThetaBin]  = Theta;
      if((sigmaT!=0.0) && (DsigmaT!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 5 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sgT_in);
    fgets(Buffer, sizeof(Buffer), File_sgT_in);
    fgets(Buffer, sizeof(Buffer), File_sgT_in);
    fgets(Buffer, sizeof(Buffer), File_sgT_in);

    sgT_pts[sgT_bin] = ThetaBin;
    sgT_en[sgT_bin] = Energy;
    sgT_bin++;
  }

  fclose(File_sgT_in);

  sprintf(Buffer, "out_%s", SvenFile);
  File_sgT_out = fopen(Buffer, "w");
  for(Int_t e=0; e<sgT_bin; e++)
  {
    fprintf(File_sgT_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", sgT_en[e], Weight, System);
    for(Int_t th=0; th<sgT_pts[e]; th++)
      fprintf(File_sgT_out, "%7.3f  %10.6f  %10.6f\n", sgT_th[e][th], sgT_val[e][th], sgT_err[e][th]);
    fprintf(File_sgT_out, "-------------------------------------------\n");
  }
  fclose(File_sgT_out);
}

//-----------------------------------------------------------------------------
