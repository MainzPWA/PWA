static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t F_val[EBINS][THBINS];
Double_t F_err[EBINS][THBINS];
Double_t F_th[EBINS][THBINS];
Double_t F_en[EBINS];
Int_t F_pts[EBINS];
Int_t F_bin;

Double_t T_val[EBINS][THBINS];
Double_t T_err[EBINS][THBINS];
Double_t T_th[EBINS][THBINS];
Double_t T_en[EBINS];
Int_t T_pts[EBINS];
Int_t T_bin;

//-----------------------------------------------------------------------------

void Parse_TF(Char_t* KashevFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBinT;
  Int_t ThetaBinF;
  Double_t Energy, Theta, CosTheta, T, DT, F, DF;
  FILE* File_TF_in;
  FILE* File_F_out;
  FILE* File_T_out;

  File_TF_in = fopen(KashevFile, "r");

  F_bin = 0;
  T_bin = 0;
  while(!feof(File_TF_in))
  {
    //Get beam energy
    if(fscanf(File_TF_in, " #       Egamma = %lf MeV\n", &Energy)!=1) break;
    ThetaBinT = 0;
    ThetaBinF = 0;
    while(fscanf(File_TF_in, "%lf %lf %lf %lf %lf %lf\n", &CosTheta, &Theta, &T, &DT, &F, &DF)==6)
    {
      F_val[F_bin][ThetaBinF] = F;
      F_err[F_bin][ThetaBinF] = DF;
      F_th[F_bin][ThetaBinF]  = Theta;
      if((F!=0.0) && (DF!=0.0)) ThetaBinF++; //Accept only 'existing' data points

      T_val[T_bin][ThetaBinT] = T;
      T_err[T_bin][ThetaBinT] = DT;
      T_th[T_bin][ThetaBinT]  = Theta;
      if((T!=0.0) && (DT!=0.0)) ThetaBinT++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_TF_in);

    F_pts[F_bin] = ThetaBinF;
    F_en[F_bin] = Energy;
    F_bin++;

    T_pts[T_bin] = ThetaBinT;
    T_en[T_bin] = Energy;
    T_bin++;
  }
  fclose(File_TF_in);

  sprintf(Buffer, "out_F_%s", KashevFile);
  File_F_out = fopen(Buffer, "w");
  for(Int_t e=0; e<F_bin; e++)
  {
    fprintf(File_F_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", F_en[e], Weight, Scale);
    for(Int_t th=0; th<F_pts[e]; th++)
      fprintf(File_F_out, "%7.3f  %10.6f  %10.6f\n", F_th[e][th], F_val[e][th], F_err[e][th]);
    fprintf(File_F_out, "-------------------------------------------\n");
  }
  fclose(File_F_out);

  sprintf(Buffer, "out_T_%s", KashevFile);
  File_T_out = fopen(Buffer, "w");
  for(Int_t e=0; e<T_bin; e++)
  {
    fprintf(File_T_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", T_en[e], Weight, System);
    for(Int_t th=0; th<T_pts[e]; th++)
      fprintf(File_F_out, "%7.3f  %10.6f  %10.6f\n", T_th[e][th], T_val[e][th], T_err[e][th]);
    fprintf(File_T_out, "-------------------------------------------\n");
  }
  fclose(File_T_out);
}

//-----------------------------------------------------------------------------
