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

void Parse_F(Char_t* PeterFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBin[512];
  Int_t Chan = 0;
  Double_t Energy, Theta, DTheta, F, DF;
  FILE* File_F_in;
  FILE* File_F_out;

  File_F_in = fopen(PeterFile, "r");
  //Skip 1 uninteresting line
  fgets(Buffer, sizeof(Buffer), File_F_in);

  for(Int_t c=0; c<512; c++) ThetaBin[c] = 0;
  while(!feof(File_F_in))
  {
    while(fscanf(File_F_in, "%d %lf %lf %lf %lf\n", &Chan, &Theta, &DTheta, &F, &DF)==5)
    {
      F_bin = Chan;
      F_val[F_bin][ThetaBin[Chan]] = F;
      F_err[F_bin][ThetaBin[Chan]] = DF;
      F_th[F_bin][ThetaBin[Chan]]  = Theta;
      F_en[F_bin] = GetBeamEnergy(Chan);
      if(DF!=0.0)
      {
        ThetaBin[Chan]++; //Accept only 'existing' data points
        F_pts[F_bin] = ThetaBin[Chan];
      }
    }
  }

  fclose(File_F_in);

  sprintf(Buffer, "out_%s", PeterFile);
  File_F_out = fopen(Buffer, "w");
  for(Int_t e=F_bin; e>-1; e--)
  {
    if(F_pts[e])
    {
      fprintf(File_F_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", F_en[e], Weight, System);
      for(Int_t th=0; th<F_pts[e]; th++)
        fprintf(File_F_out, "%7.3f  %10.6f  %10.6f\n", F_th[e][th], F_val[e][th], F_err[e][th]);
      fprintf(File_F_out, "-------------------------------------------\n");
    }
  }
  fclose(File_F_out);
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
