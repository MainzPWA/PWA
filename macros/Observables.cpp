
static const Int_t LBINS  = 10;
static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

static const Double_t MASS_PROTON   = 938.2720; //MeV
static const Double_t MASS_NEUTRON  = 939.5654; //MeV
static const Double_t MASS_PIZERO   = 134.9766; //MeV
static const Double_t MASS_PIPLUS   = 139.5702; //MeV
static const Double_t MASS_ETA      = 547.8620; //MeV
static const Double_t MASS2_PROTON  = MASS_PROTON*MASS_PROTON;
static const Double_t MASS2_NEUTRON = MASS_NEUTRON*MASS_NEUTRON;
static const Double_t MASS2_PIZERO  = MASS_PIZERO*MASS_PIZERO;
static const Double_t MASS2_PIPLUS  = MASS_PIPLUS*MASS_PIPLUS;
static const Double_t MASS2_ETA     = MASS_ETA*MASS_ETA;
static const Double_t HBARC = 197.3269; //MeV fm (hbar*c)
static const Double_t UNIT  = (HBARC/MASS_PIPLUS)*(HBARC/MASS_PIPLUS)*1e-2;

//-----------------------------------------------------------------------------

Double_t MASS_MESON  = MASS_PIZERO;
Double_t MASS2_MESON = MASS_MESON*MASS_MESON;

Int_t L_MAX;

//-----------------------------------------------------------------------------


TComplex Ep[LBINS];
TComplex Em[LBINS];
TComplex Mp[LBINS];
TComplex Mm[LBINS];

TComplex fits_Ep[LBINS][EBINS];
TComplex fits_Em[LBINS][EBINS];
TComplex fits_Mp[LBINS][EBINS];
TComplex fits_Mm[LBINS][EBINS];
Double_t fits_en[EBINS];
Int_t fits_bin;

Double_t obs_val[EBINS][THBINS];
Double_t obs_err[EBINS][THBINS];
Double_t obs_th[EBINS][THBINS];
Double_t obs_en[EBINS];
Int_t obs_pts[EBINS];
Int_t obs_bin;

//-----------------------------------------------------------------------------

Double_t L(Int_t l, Double_t x)
{
  if(l==0) return 1.0; //Legendre polynom P_0
  if(l==1) return x;   //Legendre polynom P_1
  if(l > 1) return (1.0/l) *((2.0*l - 1.0)*x*L(l-1, x) - (l - 1.0)*L(l-2, x));
}

//-----------------------------------------------------------------------------

Double_t DL(Int_t l, Double_t x)
{
  if(l==0) return 0.0; //Legendre polynom 1st derivative P'_0
  if(l==1) return 1.0; //Legendre polynom 1st derivative P'_1
  if(l > 1) return (1.0/(l - 1.0)) * ((2.0*l - 1.0)*x*DL(l-1, x) - l*DL(l-2, x));
}

//-----------------------------------------------------------------------------

Double_t D2L(Int_t l, Double_t x)
{
  if(l==0) return 0.0; //Legendre polynom 2nd derivative P"_0
  if(l==1) return 0.0; //Legendre polynom 2nd derivative P"_1
  if(l==2) return 3.0; //Legendre polynom 2nd derivative P"_2
  if(l > 2) return (1.0/(l - 2.0)) * ((2.0*l - 1.0)*x*D2L(l-1, x) - (l + 1.0)*D2L(l-2, x));
}

//-----------------------------------------------------------------------------

Double_t rho_0(Double_t omega)
{
  Double_t q = q_0(omega);
  Double_t k = omega_cm(omega);

  return (q/k);
}

//-----------------------------------------------------------------------------

Double_t q_0(Double_t omega_lab)
{
  Double_t s = MASS2_PROTON + 2.0*omega_lab*MASS_PROTON;
  Double_t sqrt_s = TMath::Sqrt(s);
  Double_t epsilon = (s + MASS2_MESON - MASS2_PROTON)/(2.0*sqrt_s);
  if(epsilon < MASS_MESON) return 0.0;
  Double_t q = TMath::Sqrt(epsilon*epsilon - MASS2_MESON);

  //printf("w_lab = %f: q_0 = %f\n", omega_lab, q);
  return q;
}

//-----------------------------------------------------------------------------

Double_t omega_cm(Double_t omega_lab)
{
  Double_t s = MASS2_PROTON + 2.0*omega_lab*MASS_PROTON;
  Double_t sqrt_s = TMath::Sqrt(s);
  Double_t omega = (s - MASS2_PROTON)/(2.0*sqrt_s);

  //printf("w_lab = %f: w_cm = %f\n", omega_lab, omega);
  return omega;
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_fit(Double_t omega)
{
  Double_t Min = 1e38;
  Int_t eM;
  for(Int_t e=0; e<fits_bin; e++)
    if(fabs(fits_en[e] - omega) < Min)
    {
      Min = fabs(fits_en[e] - omega);
      eM = e;
    }
  return eM;
}

//-----------------------------------------------------------------------------

void GetMultipoles(Double_t omega)
{
  Int_t eM = GetEnergyBin_fit(omega);

  //Pick up s,p wave multipoles from model
  Ep[0] = fits_Ep[0][eM];
  Ep[1] = fits_Ep[1][eM];
  Mp[1] = fits_Mp[1][eM];
  Mm[1] = fits_Mm[1][eM];
  //Pick up  d,f,... wave multipoles from model
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    Ep[l] = fits_Ep[l][eM];
    Em[l] = fits_Em[l][eM];
    Mp[l] = fits_Mp[l][eM];
    Mm[l] = fits_Mm[l][eM];
  }
}

//-----------------------------------------------------------------------------

Int_t GetEnergyBin_obs(Double_t omega)
{
  Double_t Min = 1e38;
  Int_t eO;

  for(Int_t e=0; e<obs_bin; e++)
    if(fabs(obs_en[e] - omega) < Min)
    {
      Min = fabs(obs_en[e] - omega);
      eO = e;
    }
  return eO;
}

//-----------------------------------------------------------------------------

TGraphErrors* GetObservable(Double_t omega)
{
  Int_t eO = GetEnergyBin_obs(omega);

  TGraphErrors* Graph = new TGraphErrors(obs_pts[eO], obs_th[eO], obs_val[eO], 0, obs_err[eO]);
  Graph->SetMarkerStyle(20);
  Graph->SetMarkerSize(0.8);

  return Graph;
}

//-----------------------------------------------------------------------------

//CGLN amplitude F1 in expansion up to l_max
TComplex F1(Double_t CosTheta)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=0; l<L_MAX+1; l++)
    Complex+=((Mp[l]*(1.0*l) + Ep[l]) * DL(l+1, CosTheta) + (Mm[l]*(1.0*l + 1.0) + Em[l]) * DL(l-1, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//CGLN amplitude F2 in expansion up to l_max
TComplex F2(Double_t CosTheta)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=1; l<L_MAX+1; l++)
    Complex+=((Mp[l]*(1.0*l + 1.0) + Mm[l]*(1.0*l)) * DL(l, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//CGLN amplitude F3 in expansion up to l_max
TComplex F3(Double_t CosTheta)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=1; l<L_MAX+1; l++)
    Complex+=((Ep[l] - Mp[l]) * D2L(l+1, CosTheta) + (Em[l] - Mm[l]) * D2L(l-1, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//CGLN amplitude F4 in expansion up to l_max
TComplex F4(Double_t CosTheta)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=2; l<L_MAX+1; l++)
    Complex+=((Mp[l] - Ep[l] - Mm[l] - Em[l]) * D2L(l, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

TComplex F1cc(Double_t CosTheta){ return TComplex::Conjugate(F1(CosTheta)); }
TComplex F2cc(Double_t CosTheta){ return TComplex::Conjugate(F2(CosTheta)); }
TComplex F3cc(Double_t CosTheta){ return TComplex::Conjugate(F3(CosTheta)); }
TComplex F4cc(Double_t CosTheta){ return TComplex::Conjugate(F4(CosTheta)); }

//-----------------------------------------------------------------------------

Double_t sigma0(Double_t theta, Double_t omega)
{
  Double_t SinTheta = TMath::Sin(theta*TMath::DegToRad());
  Double_t CosTheta = TMath::Cos(theta*TMath::DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F1(CosTheta).Rho2() + F2(CosTheta).Rho2()
                   + Sin2Theta*(0.5*F3(CosTheta).Rho2() + 0.5*F4(CosTheta).Rho2() + F2cc(CosTheta)*F3(CosTheta) + F1cc(CosTheta)*F4(CosTheta) + CosTheta*F3cc(CosTheta)*F4(CosTheta))
                   - 2.0*CosTheta*F1cc(CosTheta)*F2(CosTheta);
  return Complex.Re()*rho_0(omega)*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaS(Double_t theta, Double_t omega)
{
  Double_t SinTheta = TMath::Sin(theta*TMath::DegToRad());
  Double_t CosTheta = TMath::Cos(theta*TMath::DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = 0.5*(F3(CosTheta).Rho2() + F4(CosTheta).Rho2()) + F2cc(CosTheta)*F3(CosTheta) + F1cc(CosTheta)*F4(CosTheta) + CosTheta*F3cc(CosTheta)*F4(CosTheta);
  return -Sin2Theta*Complex.Re()*rho_0(omega)*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaT(Double_t theta, Double_t omega)
{
  Double_t SinTheta = TMath::Sin(theta*TMath::DegToRad());
  Double_t CosTheta = TMath::Cos(theta*TMath::DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F1cc(CosTheta)*F3(CosTheta) - F2cc(CosTheta)*F4(CosTheta)
                   + CosTheta*(F1cc(CosTheta)*F4(CosTheta) - F2cc(CosTheta)*F3(CosTheta))
                   - Sin2Theta*F3cc(CosTheta)*F4(CosTheta);
  return SinTheta*Complex.Im()*rho_0(omega)*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaP(Double_t theta, Double_t omega)
{
  Double_t SinTheta = TMath::Sin(theta*TMath::DegToRad());
  Double_t CosTheta = TMath::Cos(theta*TMath::DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = 2.0*F1cc(CosTheta)*F2(CosTheta) + F1cc(CosTheta)*F3(CosTheta) - F2cc(CosTheta)*F4(CosTheta)
                   - CosTheta*(F2cc(CosTheta)*F3(CosTheta) - F1cc(CosTheta)*F4(CosTheta))
                   - Sin2Theta*F3cc(CosTheta)*F4(CosTheta);
  return -SinTheta*Complex.Im()*rho_0(omega)*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaH(Double_t theta, Double_t omega)
{
  Double_t SinTheta = TMath::Sin(theta*TMath::DegToRad());
  Double_t CosTheta = TMath::Cos(theta*TMath::DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = 2.0*F1cc(CosTheta)*F2(CosTheta) + F1cc(CosTheta)*F3(CosTheta) - F2cc(CosTheta)*F4(CosTheta)
                   + CosTheta*(F1cc(CosTheta)*F4(CosTheta) - F2cc(CosTheta)*F3(CosTheta));
  return SinTheta*Complex.Im()*rho_0(omega)*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaG(Double_t theta, Double_t omega)
{
  Double_t SinTheta = TMath::Sin(theta*TMath::DegToRad());
  Double_t CosTheta = TMath::Cos(theta*TMath::DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F2cc(CosTheta)*F3(CosTheta) + F1cc(CosTheta)*F4(CosTheta);
  return Sin2Theta*Complex.Im()*rho_0(omega)*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaF(Double_t theta, Double_t omega)
{
  Double_t SinTheta = TMath::Sin(theta*TMath::DegToRad());
  Double_t CosTheta = TMath::Cos(theta*TMath::DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F1cc(CosTheta)*F3(CosTheta) - F2cc(CosTheta)*F4(CosTheta)
                   - CosTheta*(F2cc(CosTheta)*F3(CosTheta) - F1cc(CosTheta)*F4(CosTheta));
  return SinTheta*Complex.Re()*rho_0(omega)*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaE(Double_t theta, Double_t omega)
{
  Double_t SinTheta = TMath::Sin(theta*TMath::DegToRad());
  Double_t CosTheta = TMath::Cos(theta*TMath::DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F1cc(CosTheta)*F1(CosTheta) + F2cc(CosTheta)*F2(CosTheta)
                   - 2.0*CosTheta*F1cc(CosTheta)*F2(CosTheta)
                   + Sin2Theta*(F2cc(CosTheta)*F3(CosTheta) + F1cc(CosTheta)*F4(CosTheta));
  return Complex.Re()*rho_0(omega)*UNIT;
}

//-----------------------------------------------------------------------------

Double_t S(Double_t theta, Double_t omega)
{
  return sigmaS(Theta, Omega)/sigma0(Theta, Omega);
}

//-----------------------------------------------------------------------------

Double_t T(Double_t theta, Double_t omega)
{
  return sigmaT(Theta, Omega)/sigma0(Theta, Omega);
}

//-----------------------------------------------------------------------------

Double_t P(Double_t theta, Double_t omega)
{
  return sigmaP(Theta, Omega)/sigma0(Theta, Omega);
}

//-----------------------------------------------------------------------------

Double_t H(Double_t theta, Double_t omega)
{
  return sigmaH(Theta, Omega)/sigma0(Theta, Omega);
}

//-----------------------------------------------------------------------------

Double_t G(Double_t theta, Double_t omega)
{
  return sigmaG(Theta, Omega)/sigma0(Theta, Omega);
}

//-----------------------------------------------------------------------------

Double_t F(Double_t theta, Double_t omega)
{
  return sigmaF(Theta, Omega)/sigma0(Theta, Omega);
}

//-----------------------------------------------------------------------------

Double_t E(Double_t theta, Double_t omega)
{
  return sigmaE(Theta, Omega)/sigma0(Theta, Omega);
}

//-----------------------------------------------------------------------------

void Parse_Fits()
{
  Char_t Buffer[1024];
  Double_t Energy, Re, DRe, Im, DIm;
  FILE* Fits_Ep[LBINS];
  FILE* Fits_Em[LBINS];
  FILE* Fits_Mp[LBINS];
  FILE* Fits_Mm[LBINS];

  //Open s,p waves
  Fits_Ep[0] = fopen("plots/E0p.txt", "r");
  Fits_Ep[1] = fopen("plots/E1p.txt", "r");
  Fits_Mp[1] = fopen("plots/M1p.txt", "r");
  Fits_Mm[1] = fopen("plots/M1m.txt", "r");
  //Open d,f,... waves
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    sprintf(Buffer, "plots/E%dp.txt", l); Fits_Ep[l] = fopen(Buffer, "r");
    sprintf(Buffer, "plots/E%dm.txt", l); Fits_Em[l] = fopen(Buffer, "r");
    sprintf(Buffer, "plots/M%dp.txt", l); Fits_Mp[l] = fopen(Buffer, "r");
    sprintf(Buffer, "plots/M%dm.txt", l); Fits_Mm[l] = fopen(Buffer, "r");
  }

  //Load E0+ multipole
  fits_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), Fits_Ep[0]);
  fgets(Buffer, sizeof(Buffer), Fits_Ep[0]);
  while(!feof(Fits_Ep[0]))
  {
    if(fscanf(Fits_Ep[0], "%lf %lf %lf %lf %lf\n", &Energy, &Re, &DRe, &Im, &DIm)!=5) break;
    fits_Ep[0][fits_bin] = TComplex(Re, Im);
    fits_en[fits_bin] = (Energy*Energy - MASS2_PROTON)/(2.0*MASS_PROTON);
    fits_bin++;
  }

  //Load E1+ multipole
  fits_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), Fits_Ep[1]);
  fgets(Buffer, sizeof(Buffer), Fits_Ep[1]);
  while(!feof(Fits_Ep[1]))
  {
    if(fscanf(Fits_Ep[1], "%lf %lf %lf %lf %lf\n", &Energy, &Re, &DRe, &Im, &DIm)!=5) break;
    fits_Ep[1][fits_bin] = TComplex(Re, Im);
    fits_en[fits_bin] = (Energy*Energy - MASS2_PROTON)/(2.0*MASS_PROTON);
    fits_bin++;
  }

  //Load M1+ multipole
  fits_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), Fits_Mp[1]);
  fgets(Buffer, sizeof(Buffer), Fits_Mp[1]);
  while(!feof(Fits_Mp[1]))
  {
    if(fscanf(Fits_Mp[1], "%lf %lf %lf %lf %lf\n", &Energy, &Re, &DRe, &Im, &DIm)!=5) break;
    fits_Mp[1][fits_bin] = TComplex(Re, Im);
    fits_en[fits_bin] = (Energy*Energy - MASS2_PROTON)/(2.0*MASS_PROTON);
    fits_bin++;
  }

  //Load M1- multipole
  fits_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), Fits_Mm[1]);
  fgets(Buffer, sizeof(Buffer), Fits_Mm[1]);
  while(!feof(Fits_Mm[1]))
  {
    if(fscanf(Fits_Mm[1], "%lf %lf %lf %lf %lf\n", &Energy, &Re, &DRe, &Im, &DIm)!=5) break;
    fits_Mm[1][fits_bin] = TComplex(Re, Im);
    fits_en[fits_bin] = (Energy*Energy - MASS2_PROTON)/(2.0*MASS_PROTON);
    fits_bin++;
  }

   //Load El+ multipole for d,f,... waves
  for(Int_t l=2;l<L_MAX+1; l++)
  {
    fits_bin = 0;
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), Fits_Ep[l]);
    fgets(Buffer, sizeof(Buffer), Fits_Ep[l]);
    while(!feof(Fits_Ep[l]))
    {
      if(fscanf(Fits_Ep[l], "%lf %lf %lf %lf %lf\n", &Energy, &Re, &DRe, &Im, &DIm)!=5) break;
      fits_Ep[l][fits_bin] = TComplex(Re, Im);
      fits_en[fits_bin] = (Energy*Energy - MASS2_PROTON)/(2.0*MASS_PROTON);
      fits_bin++;
    }
  }

   //Load El- multipole for d,f,... waves
  for(Int_t l=2;l<L_MAX+1; l++)
  {
    fits_bin = 0;
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), Fits_Em[l]);
    fgets(Buffer, sizeof(Buffer), Fits_Em[l]);
    while(!feof(Fits_Em[l]))
    {
      if(fscanf(Fits_Em[l], "%lf %lf %lf %lf %lf\n", &Energy, &Re, &DRe, &Im, &DIm)!=5) break;
      fits_Em[l][fits_bin] = TComplex(Re, Im);
      fits_en[fits_bin] = (Energy*Energy - MASS2_PROTON)/(2.0*MASS_PROTON);
      fits_bin++;
    }
  }

   //Load Ml+ multipole for d,f,... waves
  for(Int_t l=2;l<L_MAX+1; l++)
  {
    fits_bin = 0;
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), Fits_Mp[l]);
    fgets(Buffer, sizeof(Buffer), Fits_Mp[l]);
    while(!feof(Fits_Mp[l]))
    {
      if(fscanf(Fits_Mp[l], "%lf %lf %lf %lf %lf\n", &Energy, &Re, &DRe, &Im, &DIm)!=5) break;
      fits_Mp[l][fits_bin] = TComplex(Re, Im);
      fits_en[fits_bin] = (Energy*Energy - MASS2_PROTON)/(2.0*MASS_PROTON);
      fits_bin++;
    }
  }

   //Load Ml- multipole for d,f,... waves
  for(Int_t l=2;l<L_MAX+1; l++)
  {
    fits_bin = 0;
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), Fits_Mm[l]);
    fgets(Buffer, sizeof(Buffer), Fits_Mm[l]);
    while(!feof(Fits_Mm[l]))
    {
      if(fscanf(Fits_Mm[l], "%lf %lf %lf %lf %lf\n", &Energy, &Re, &DRe, &Im, &DIm)!=5) break;
      fits_Mm[l][fits_bin] = TComplex(Re, Im);
      fits_en[fits_bin] = (Energy*Energy - MASS2_PROTON)/(2.0*MASS_PROTON);
      fits_bin++;
    }
  }

  //Close s,p waves
  fclose(Fits_Ep[0]);
  fclose(Fits_Ep[1]);
  fclose(Fits_Mp[1]);
  fclose(Fits_Mm[1]);
  //Close d,f,... waves
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    fclose(Fits_Ep[l]);
    fclose(Fits_Em[l]);
    fclose(Fits_Mp[l]);
    fclose(Fits_Mm[l]);
  }

  return;

  //Debug output
  printf("EBins: %d\n", fits_bin);
  for(Int_t e=0; e<fits_bin; e++)
  {
    printf("%d (%f MeV)\n", e, fits_en[e]);
    printf("%f %f %f %f %f %f %f %f\n",
           fits_Ep[0][e].Re(), fits_Ep[1][e].Re(), fits_Mp[1][e].Re(), fits_Mm[1][e].Re(), fits_Ep[2][e].Re(), fits_Mp[2][e].Re(), fits_Em[2][e].Re(), fits_Mm[2][e].Re());
    printf("%f %f %f %f %f %f %f %f\n",
           fits_Ep[0][e].Im(), fits_Ep[1][e].Im(), fits_Mp[1][e].Im(), fits_Mm[1][e].Im(), fits_Ep[2][e].Im(), fits_Mp[2][e].Im(), fits_Em[2][e].Im(), fits_Mm[2][e].Im());
  }
}

//-----------------------------------------------------------------------------

void Parse_Obsv(Char_t* Filename)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Double_t Energy, Weight, Scale, Theta, Obs, DObs;
  FILE* File_obs;

  File_obs = fopen(Filename, "r");

  obs_bin = 0;
  while(!feof(File_obs))
  {
    //Get beam energy
    if(fscanf(File_obs, "E = %lf MeV, Wt = %lf, Sc = %lf\n", &Energy, &Weight, &Scale)!=3) break;
    ThetaBin = 0;
    while(fscanf(File_obs, "%lf %lf %lf\n", &Theta, &Obs, &DObs)==3)
    {
      obs_val[obs_bin][ThetaBin] = Obs;
      obs_err[obs_bin][ThetaBin] = DObs;
      obs_th[obs_bin][ThetaBin]  = Theta;
      if((Obs!=0.0) && (DObs!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_obs);

    obs_pts[obs_bin] = ThetaBin;
    obs_en[obs_bin] = Energy;
    obs_bin++;
  }

  fclose(File_obs);
  return;

  //Debug output
  printf("EBins: %d\n", obs_bin);
  for(Int_t e=0; e<obs_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, obs_en[e], obs_pts[e]);
    for(Int_t th=0; th<obs_pts[e]; th++)
      printf("%f %f %f\n", obs_th[e][th], obs_val[e][th], obs_err[e][th]);
  }
}

//-----------------------------------------------------------------------------

TH1F* Draw_sigma0(Double_t omega, Int_t lmax = 3)
{
  L_MAX = lmax;
  Parse_Fits();
  GetMultipoles(omega);
  Parse_Obsv("data/sg0.txt");
  TGraphErrors* Graph = GetObservable(omega);

  TH1F* Histo = new TH1F("sigma0", "sigma0", 181, -0.5, 180.5);
  for(Int_t th=0; th<181; th++)
    Histo->SetBinContent(th+1, sigma0(th*1.0, omega));
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#sigma_{0} / #mub/sr");

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_sigmaS(Double_t omega, Int_t lmax = 3)
{
  L_MAX = lmax;
  Parse_Fits();
  GetMultipoles(omega);
  Parse_Obsv("data/sgS.txt");
  TGraphErrors* Graph = GetObservable(omega);

  TH1F* Histo = new TH1F("sigmaS", "sigmaS", 181, -0.5, 180.5);
  for(Int_t th=0; th<181; th++)
    Histo->SetBinContent(th+1, sigmaS(th*1.0, omega));
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#sigma_{#Sigma} / #mub/sr");

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_sigmaT(Double_t omega, Int_t lmax = 3)
{
  L_MAX = lmax;
  Parse_Fits();
  GetMultipoles(omega);
  Parse_Obsv("data/sgT.txt");
  TGraphErrors* Graph = GetObservable(omega);

  TH1F* Histo = new TH1F("sigmaT", "sigmaT", 181, -0.5, 180.5);
  for(Int_t th=0; th<181; th++)
    Histo->SetBinContent(th+1, sigmaT(th*1.0, omega));
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#sigma_{T} / #mub/sr");

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_sigmaP(Double_t omega, Int_t lmax = 3)
{
  L_MAX = lmax;
  Parse_Fits();
  GetMultipoles(omega);
  Parse_Obsv("data/sgP.txt");
  TGraphErrors* Graph = GetObservable(omega);

  TH1F* Histo = new TH1F("sigmaP", "sigmaP", 181, -0.5, 180.5);
  for(Int_t th=0; th<181; th++)
    Histo->SetBinContent(th+1, sigmaP(th*1.0, omega));
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#sigma_{P} / #mub/sr");

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_sigmaE(Double_t omega, Int_t lmax = 3)
{
  L_MAX = lmax;
  Parse_Fits();
  GetMultipoles(omega);
  Parse_Obsv("data/sgE.txt");
  TGraphErrors* Graph = GetObservable(omega);

  TH1F* Histo = new TH1F("sigmaE", "sigmaE", 181, -0.5, 180.5);
  for(Int_t th=0; th<181; th++)
    Histo->SetBinContent(th+1, sigmaE(th*1.0, omega));
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#sigma_{E} / #mub/sr");

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_sigmaF(Double_t omega, Int_t lmax = 3)
{
  L_MAX = lmax;
  Parse_Fits();
  GetMultipoles(omega);
  Parse_Obsv("data/sgF.txt");
  TGraphErrors* Graph = GetObservable(omega);

  TH1F* Histo = new TH1F("sigmaF", "sigmaF", 181, -0.5, 180.5);
  for(Int_t th=0; th<181; th++)
    Histo->SetBinContent(th+1, sigmaF(th*1.0, omega));
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#sigma_{F} / #mub/sr");

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_sigmaG(Double_t omega, Int_t lmax = 3)
{
  L_MAX = lmax;
  Parse_Fits();
  GetMultipoles(omega);
  Parse_Obsv("data/sgG.txt");
  TGraphErrors* Graph = GetObservable(omega);

  TH1F* Histo = new TH1F("sigmaG", "sigmaG", 181, -0.5, 180.5);
  for(Int_t th=0; th<181; th++)
    Histo->SetBinContent(th+1, sigmaG(th*1.0, omega));
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#sigma_{G} / #mub/sr");

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_sigmaH(Double_t omega, Int_t lmax = 3)
{
  L_MAX = lmax;
  Parse_Fits();
  GetMultipoles(omega);
  Parse_Obsv("data/sgH.txt");
  TGraphErrors* Graph = GetObservable(omega);

  TH1F* Histo = new TH1F("sigmaH", "sigmaH", 181, -0.5, 180.5);
  for(Int_t th=0; th<181; th++)
    Histo->SetBinContent(th+1, sigmaH(th*1.0, omega));
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#sigma_{H} / #mub/sr");

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_S(Double_t omega, Int_t lmax = 3)
{
  TH1F* HistoS = Draw_sigmaS(omega, lmax);
  TH1F* Histo0 = Draw_sigma0(omega, lmax);

  TH1F* Histo = (TH1F*)HistoS->Clone("S");
  Histo->Divide(Histo0);
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("#Sigma");
  Histo->SetMinimum(-1.0);
  Histo->SetMaximum(+1.0);

  Parse_Obsv("data/S.txt");
  TGraphErrors* Graph = GetObservable(omega);

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_T(Double_t omega, Int_t lmax = 3)
{
  TH1F* HistoT = Draw_sigmaT(omega, lmax);
  TH1F* Histo0 = Draw_sigma0(omega, lmax);

  TH1F* Histo = (TH1F*)HistoT->Clone("T");
  Histo->Divide(Histo0);
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("T");
  Histo->SetMinimum(-1.0);
  Histo->SetMaximum(+1.0);

  Parse_Obsv("data/T.txt");
  TGraphErrors* Graph = GetObservable(omega);

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_P(Double_t omega, Int_t lmax = 3)
{
  TH1F* HistoP = Draw_sigmaP(omega, lmax);
  TH1F* Histo0 = Draw_sigma0(omega, lmax);

  TH1F* Histo = (TH1F*)HistoP->Clone("P");
  Histo->Divide(Histo0);
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("P");
  Histo->SetMinimum(-1.0);
  Histo->SetMaximum(+1.0);

  Parse_Obsv("data/P.txt");
  TGraphErrors* Graph = GetObservable(omega);

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_E(Double_t omega, Int_t lmax = 3)
{
  TH1F* HistoE = Draw_sigmaE(omega, lmax);
  TH1F* Histo0 = Draw_sigma0(omega, lmax);

  TH1F* Histo = (TH1F*)HistoE->Clone("E");
  Histo->Divide(Histo0);
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("E");
  Histo->SetMinimum(-1.0);
  Histo->SetMaximum(+1.0);

  Parse_Obsv("data/E.txt");
  TGraphErrors* Graph = GetObservable(omega);

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_F(Double_t omega, Int_t lmax = 3)
{
  TH1F* HistoF = Draw_sigmaF(omega, lmax);
  TH1F* Histo0 = Draw_sigma0(omega, lmax);

  TH1F* Histo = (TH1F*)HistoF->Clone("F");
  Histo->Divide(Histo0);
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("F");
  Histo->SetMinimum(-1.0);
  Histo->SetMaximum(+1.0);

  Parse_Obsv("data/F.txt");
  TGraphErrors* Graph = GetObservable(omega);

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_G(Double_t omega, Int_t lmax = 3)
{
  TH1F* HistoG = Draw_sigmaG(omega, lmax);
  TH1F* Histo0 = Draw_sigma0(omega, lmax);

  TH1F* Histo = (TH1F*)HistoG->Clone("G");
  Histo->Divide(Histo0);
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("G");
  Histo->SetMinimum(-1.0);
  Histo->SetMaximum(+1.0);

  Parse_Obsv("data/G.txt");
  TGraphErrors* Graph = GetObservable(omega);

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------

TH1F* Draw_H(Double_t omega, Int_t lmax = 3)
{
  TH1F* HistoH = Draw_sigmaH(omega, lmax);
  TH1F* Histo0 = Draw_sigma0(omega, lmax);

  TH1F* Histo = (TH1F*)HistoH->Clone("H");
  Histo->Divide(Histo0);
  Histo->GetXaxis()->SetTitle("#theta / deg");
  Histo->GetYaxis()->SetTitle("H");
  Histo->SetMinimum(-1.0);
  Histo->SetMaximum(+1.0);

  Parse_Obsv("data/H.txt");
  TGraphErrors* Graph = GetObservable(omega);

  Histo->Draw("c");
  Graph->Draw("pz");

  return Histo;
}

//-----------------------------------------------------------------------------
