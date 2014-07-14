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

void Parse_T(Char_t* PeterFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBin[512];
  Int_t Chan = 0;
  Double_t Energy, Theta, DTheta, T, DT;
  FILE* File_T_in;
  FILE* File_T_out;

  File_T_in = fopen(PeterFile, "r");
  //Skip 1 uninteresting line
  fgets(Buffer, sizeof(Buffer), File_T_in);

  for(Int_t c=0; c<512; c++) ThetaBin[c] = 0;
  while(!feof(File_T_in))
  {
    while(fscanf(File_T_in, "%d %lf %lf %lf %lf\n", &Chan, &Theta, &DTheta, &T, &DT)==5)
    {
      T_bin = Chan;
      T_val[T_bin][ThetaBin[Chan]] = T;
      T_err[T_bin][ThetaBin[Chan]] = DT;
      T_th[T_bin][ThetaBin[Chan]]  = Theta;
      T_en[T_bin] = GetBeamEnergy(Chan);
      if((T!=0.0) && (DT!=0.0))
      {
        ThetaBin[Chan]++; //Accept only 'existing' data points
        T_pts[T_bin] = ThetaBin[Chan];
      }
    }
  }

  fclose(File_T_in);

  sprintf(Buffer, "out_%s", PeterFile);
  File_T_out = fopen(Buffer, "w");
  for(Int_t e=T_bin; e>-1; e--)
  {
    if(T_pts[e])
    {
      fprintf(File_T_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", T_en[e], Weight, System);
      for(Int_t th=0; th<T_pts[e]; th++)
        fprintf(File_T_out, "%7.3f  %10.6f  %10.6f\n", T_th[e][th], T_val[e][th], T_err[e][th]);
      fprintf(File_T_out, "-------------------------------------------\n");
    }
  }
  fclose(File_T_out);
}

//-----------------------------------------------------------------------------

Double_t GetBeamEnergy(Int_t chan, Double_t* wdth=NULL)
{
  Double_t ChanEn[352];
  Double_t ChanWd[352];

  FILE* ChanEnFile = fopen("ChanEn450.txt", "r");
  for(Int_t t=0; t<352; t++)
    fscanf(ChanEnFile, "%lf %lf", &ChanEn[t], &ChanWd[t]);
  fclose(ChanEnFile);

  Double_t ener = ChanEn[chan];
  Double_t widt = ChanWd[chan]/2.0;

  //printf("Channel %3d = %6.2f +- %4.2f MeV\n", chan, ener, widt);
  if(wdth) (*wdth) = widt;
  return ener;
}

//-----------------------------------------------------------------------------
