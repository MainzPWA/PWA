static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t T_val[EBINS][THBINS];
Double_t T_err[EBINS][THBINS];
Double_t T_th[EBINS][THBINS];
Double_t T_en[EBINS];
Int_t T_pts[EBINS];
Int_t T_bin;

//-----------------------------------------------------------------------------

void Parse_T(Char_t* KashevFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Double_t Energy, CosTheta, T, DT;
  FILE* File_T_in;
  FILE* File_T_out;

  File_T_in = fopen(KashevFile, "r");

  T_bin = 0;
  while(!feof(File_T_in))
  {
    //Get beam energy
    if(fscanf(File_T_in, "#       E_gamma = %lf MeV\n", &Energy)!=1) break;
    ThetaBin = 0;
    while(fscanf(File_T_in, "%lf %lf %lf\n", &CosTheta, &T, &DT)==3)
    {
      T_val[T_bin][ThetaBin] = T;
      T_err[T_bin][ThetaBin] = DT;
      T_th[T_bin][ThetaBin]  = TMath::ACos(CosTheta)*TMath::RadToDeg();
      if((T!=0.0) && (DT!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_T_in);

    T_pts[T_bin] = ThetaBin;
    T_en[T_bin] = Energy;
    T_bin++;
  }

  fclose(File_T_in);

  sprintf(Buffer, "out_%s", KashevFile);
  File_T_out = fopen(Buffer, "w");
  for(Int_t e=0; e<T_bin; e++)
  {
    fprintf(File_T_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", T_en[e], Weight, System);
    for(Int_t th=0; th<T_pts[e]; th++)
      fprintf(File_T_out, "%7.3f  %10.6f  %10.6f\n", T_th[e][th], T_val[e][th], T_err[e][th]);
    fprintf(File_T_out, "-------------------------------------------\n");
  }
  fclose(File_T_out);
}

//-----------------------------------------------------------------------------
