static const Int_t EBINS  = 1024;
static const Int_t THBINS = 128;
static const Double_t D2R = TMath::DegToRad();

Double_t sg0_val[EBINS][THBINS];
Double_t sg0_err[EBINS][THBINS];
Double_t sg0_th[EBINS][THBINS];
Double_t sg0_en[EBINS];
Double_t sg0_wt[EBINS];
Double_t sg0_sy[EBINS];
Int_t sg0_pts[EBINS];
Int_t sg0_bin;

//-----------------------------------------------------------------------------

void Parse_sg0()
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Double_t Energy, Weight, System, Theta, sigma0, Dsigma0;
  FILE* File_sg0;

  File_sg0 = fopen("sg0.txt", "r");

  sg0_bin = 0;
  while(!feof(File_sg0))
  {
    //Get beam energy
    if(fscanf(File_sg0, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;
    ThetaBin = 0;
    while(fscanf(File_sg0, "%lf %lf %lf\n", &Theta, &sigma0, &Dsigma0)==3)
    {
      sg0_val[sg0_bin][ThetaBin] = sigma0;
      sg0_err[sg0_bin][ThetaBin] = Dsigma0;
      sg0_th[sg0_bin][ThetaBin]  = Theta;
      ThetaBin++;
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sg0);

    sg0_pts[sg0_bin] = ThetaBin;
    sg0_en[sg0_bin] = Energy;
    sg0_wt[sg0_bin] = Weight;
    sg0_sy[sg0_bin] = System;
    sg0_bin++;
  }

  fclose(File_sg0);
  return;

  //Debug output
  printf("EBins: %d\n", sg0_bin);
  for(Int_t e=0; e<sg0_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, sg0_en[e], sg0_pts[e]);
    for(Int_t th=0; th<sg0_pts[e]; th++)
      printf("%f %f %f\n", sg0_th[e][th], sg0_val[e][th], sg0_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

TH1F* sg0(Double_t Energy0, Int_t Bins=18)
{
  Int_t Bin0;
  Char_t Buffer[256];
  TH1F* sg0Data;
  TH1F* sg0Fit;
  TF1* FitABC;
  TFitResult* FitResult;
  Double_t BinWidth0;
  Int_t BinNumber;
  Double_t ThetaLo, ThetaHi;
  Double_t CosThetaLo, CosThetaHi;
  Double_t sigma, Dsigma;
  Double_t Param[3];
  Double_t* Error;
  Double_t Chi2;
  Int_t NDF;

  //Load data
  Parse_sg0();
  Bin0 = GetEnergyBin_sg0(Energy0);
  Energy0 = sg0_en[Bin0];

  //Create sg0 histogram in cos(theta) domain
  sprintf(Buffer, "Data_sg0_%3.0fMeV", Energy0);
  sg0Data = new TH1F(Buffer, Buffer, sg0_pts[Bin0], -1.0, 1.0);
  BinWidth0 = 2.0/sg0_pts[Bin0];
  for(Int_t n=0; n<sg0_pts[Bin0]; n++)
  //Fill sg0 values and errors
  {
    BinNumber = 1 + (cos(sg0_th[Bin0][n]*D2R)+1.0)/BinWidth0;
    sg0Data->SetBinContent(BinNumber, sg0_val[Bin0][n]);
    sg0Data->SetBinError(BinNumber, sg0_err[Bin0][n]);
  }

  //Fit ABC coefficients to cross section data
  FitABC = new TF1("ABC", ABC, -1.0, 1.0, 3);
  TFitResult* FitResult = (sg0Data->Fit(FitABC, "RQS")).Get();  //Perform fit
  FitABC->GetParameters(Param);   //Get back fitted parameters
  Error = FitABC->GetParErrors(); //Get back errors on fit parameters
  Chi2 = FitABC->GetChisquare();
  NDF = FitABC->GetNDF();
  //printf("A = %f, B = %f, C = %f\n", Param[0], Param[1], Param[2]);

  //Create fitted sg0 histogram in theta domain
  sprintf(Buffer, "Fit_sg0_%3.0fMeV", Energy0);
  sg0Fit = new TH1F(Buffer, Buffer, Bins, 0.0, 180.0);
  BinWidth0 = 180.0/Bins;
  //Fill sg0 values and errors calculated from fit function
  for(Int_t n=0; n<Bins; n++)
  {
    BinNumber = n+1;
    ThetaLo = BinWidth0*n;
    ThetaHi = BinWidth0*n + BinWidth0;
    CosThetaLo = cos(ThetaLo*D2R);
    CosThetaHi = cos(ThetaHi*D2R);
    //The following lines use the Integral() and IntegralError() methods to calculate the
    //function value over the bin range with correct errors, including the covariance matrix(!)
    sigma = FitABC->Integral(CosThetaLo, CosThetaHi)/(CosThetaHi-CosThetaLo);
    Dsigma = FitABC->IntegralError(CosThetaLo, CosThetaHi, FitResult->GetParams(), FitResult->GetCovarianceMatrix().GetMatrixArray())/(CosThetaHi-CosThetaLo);
    sg0Fit->SetBinContent(BinNumber, sigma);
    sg0Fit->SetBinError(BinNumber, Dsigma);
  }

  //Clean up intermediate objects
  sg0Data->Delete();
  FitABC->Delete();

  return sg0Fit;
}

//-----------------------------------------------------------------------------

TH1F* sg0_cos(Double_t Energy0, Int_t Bins=20)
{
  Int_t Bin0;
  Char_t Buffer[256];
  TH1F* sg0Data;
  TH1F* sg0Fit;
  TF1* FitABC;
  TFitResult* FitResult;
  Double_t BinWidth0;
  Int_t BinNumber;
  Double_t ThetaLo, ThetaHi;
  Double_t CosThetaLo, CosThetaHi;
  Double_t sigma, Dsigma;
  Double_t Param[3];
  Double_t* Error;
  Double_t Chi2;
  Int_t NDF;

  //Load data
  Parse_sg0();
  Bin0 = GetEnergyBin_sg0(Energy0);
  Energy0 = sg0_en[Bin0];

  //Create sg0 histogram in cos(theta) domain
  sprintf(Buffer, "Data_sg0_%3.0fMeV", Energy0);
  sg0Data = new TH1F(Buffer, Buffer, sg0_pts[Bin0], -1.0, 1.0);
  BinWidth0 = 2.0/sg0_pts[Bin0];
  for(Int_t n=0; n<sg0_pts[Bin0]; n++)
  //Fill sg0 values and errors
  {
    BinNumber = 1 + (cos(sg0_th[Bin0][n]*D2R)+1.0)/BinWidth0;
    sg0Data->SetBinContent(BinNumber, sg0_val[Bin0][n]);
    sg0Data->SetBinError(BinNumber, sg0_err[Bin0][n]);
  }

  //Fit ABC coefficients to cross section data
  FitABC = new TF1("ABC", ABC, -1.0, 1.0, 3);
  TFitResult* FitResult = (sg0Data->Fit(FitABC, "RQS")).Get();  //Perform fit
  FitABC->GetParameters(Param);   //Get back fitted parameters
  Error = FitABC->GetParErrors(); //Get back errors on fit parameters
  Chi2 = FitABC->GetChisquare();
  NDF = FitABC->GetNDF();
  //printf("A = %f, B = %f, C = %f\n", Param[0], Param[1], Param[2]);

  //Create fitted sg0 histogram in cos(theta) domain
  sprintf(Buffer, "Fit_sg0_%3.0fMeV", Energy0);
  sg0Fit = new TH1F(Buffer, Buffer, Bins, -1.0, 1.0);
  BinWidth0 = 2.0/Bins;
  //Fill sg0 values and errors calculated from fit function
  for(Int_t n=0; n<Bins; n++)
  {
    BinNumber = n+1;
    CosThetaLo = (BinWidth0*n) - 1.0;
    CosThetaHi = (BinWidth0*n + BinWidth0) - 1.0;
    //The following lines use the Integral() and IntegralError() methods to calculate the
    //function value over the bin range with correct errors, including the covariance matrix(!)
    sigma = FitABC->Integral(CosThetaLo, CosThetaHi)/(CosThetaHi-CosThetaLo);
    Dsigma = FitABC->IntegralError(CosThetaLo, CosThetaHi, FitResult->GetParams(), FitResult->GetCovarianceMatrix().GetMatrixArray())/(CosThetaHi-CosThetaLo);
    sg0Fit->SetBinContent(BinNumber, sigma);
    sg0Fit->SetBinError(BinNumber, Dsigma);
  }

  //Clean up intermediate objects
  sg0Data->Delete();
  FitABC->Delete();

  return sg0Fit;
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_sg0(Double_t pEnergy)
{
  //Get energy bin for sigma0 for given energy
  Double_t Min = 1e38;
  Int_t e0;

  for(Int_t e=0; e<sg0_bin; e++)
    if(fabs(sg0_en[e] - pEnergy) < Min)
    {
      Min = fabs(sg0_en[e] - pEnergy);
      e0 = e;
    }
  return e0;
}

//-----------------------------------------------------------------------------

//Define fit function A + B*x + C*x^2
Double_t ABC(Double_t* x, Double_t* Param)
{
  //Param[0] = A
  //Param[1] = B
  //Param[2] = B
  //x[0] = running variable (=cos theta)
  return (Param[0] + Param[1]*x[0] + Param[2]*x[0]*x[0]);
}

//-----------------------------------------------------------------------------

void Interpolate_sg0(Int_t Bins=18, Double_t Weight=1.00, Double_t System=0.05)
{
  FILE* Source;
  FILE* Target;
  FILE* Output;
  Double_t omega_Source[EBINS];
  Double_t omega_Target[EBINS];
  Double_t Value, Error, Width;
  Int_t N_Source = 0;
  Int_t N_Target = 0;
  TH1F* sg0Fit;

  //Load source energy list
  Source = fopen("omega_source.txt", "r");
  while(!feof(Source))
  {
    if(fscanf(Source, "%lf", &Value)!=1) break;
    omega_Source[N_Source] = Value;
    N_Source++;
  }
  fclose(Source);

  //Load target energy list
  Target = fopen("omega_target.txt", "r");
  while(!feof(Target))
  {
    if(fscanf(Target, "%lf", &Value)!=1) break;
    omega_Target[N_Target] = Value;
    N_Target++;
  }
  fclose(Target);

  //Create interpolators for values and errors
  TMultiDim fint_sg0(2);
  TMultiDim fint_Dsg0(2);
  fint_sg0.SetNBins(N_Source, Bins);
  fint_Dsg0.SetNBins(N_Source, Bins);

  //Fill interpolators with cross section values and errors at source energies
  Width = 180.0/Bins;
  for(Int_t n=0; n<N_Source; n++)
  {
    sg0Fit = sg0(omega_Source[n], Bins);
    fint_sg0.SetCoord(0, n, omega_Source[n]);
    fint_Dsg0.SetCoord(0, n, omega_Source[n]);
    for(Int_t t=0; t<Bins; t++)
    {
      fint_sg0.SetValue(n, t, sg0Fit->GetBinContent(t+1));
      fint_sg0.SetCoord(1, t, Width*t + 0.5*Width);
      fint_Dsg0.SetValue(n, t, sg0Fit->GetBinError(t+1)*sg0Fit->GetBinError(t+1));
      fint_Dsg0.SetCoord(1, t, Width*t + 0.5*Width);
    }
    sg0Fit->Delete();
  }

  //Get cross section values and errors at target energies
  Output = fopen("sg0_interpolated.txt", "w");
  for(Int_t n=0; n<N_Target; n++)
  {
    fprintf(Output, "E = %8.3f MeV, Wght = %4.2f, Syst = %4.2f\n", omega_Target[n], Weight, System);
    for(Int_t t=0; t<Bins; t++)
    {
      Value = fint_sg0.Interpolate(omega_Target[n], Width*t + 0.5*Width);
      Error = TMath::Sqrt(fabs(fint_Dsg0.Interpolate(omega_Target[n], Width*t + 0.5*Width)));
      fprintf(Output, "%7.3f  %10.6f  %10.6f\n", Width*t + 0.5*Width, Value, Error);
    }
    fprintf(Output, "------------------------------------------\n");
  }
  fclose(Output);
}

//-----------------------------------------------------------------------------

