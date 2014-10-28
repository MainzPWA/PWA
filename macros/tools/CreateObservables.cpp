using namespace TMath;

//-----------------------------------------------------------------------------

static const Int_t LBINS = 11;
static const Int_t EBINS = 1024;
static const Int_t L_MAX = 10;
static const Double_t MASS_MESON    =  493.6770; //MeV
static const Double_t MASS_INITIAL  =  938.2720; //MeV
static const Double_t MASS_FINAL    = 1115.6830; //MeV
static const Double_t MASS2_MESON   = MASS_MESON*MASS_MESON;
static const Double_t MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;
static const Double_t MASS2_FINAL   = MASS_FINAL*MASS_FINAL;
static const Double_t MASS_PIPLUS   = 139.5702; //MeV
static const Double_t HBARC         = 197.3269; //MeV fm (hbar*c)
static const Double_t UNIT          = (HBARC/MASS_PIPLUS)*(HBARC/MASS_PIPLUS)*1e-2;

//-----------------------------------------------------------------------------

TComplex maid_Ep[LBINS][EBINS];
TComplex maid_Em[LBINS][EBINS];
TComplex maid_Mp[LBINS][EBINS];
TComplex maid_Mm[LBINS][EBINS];
Double_t maid_en[EBINS];
Int_t maid_bin;

//-----------------------------------------------------------------------------

Double_t q_mes(Double_t omega_lab)
{
  Double_t s = MASS2_INITIAL + 2.0*omega_lab*MASS_INITIAL;
  Double_t sqrt_s = Sqrt(s);
  Double_t epsilon = (s + MASS2_MESON - MASS2_FINAL)/(2.0*sqrt_s);
  if(epsilon < MASS_MESON) return 0.0;
  Double_t q = Sqrt(epsilon*epsilon - MASS2_MESON);

  return q;
}

//-----------------------------------------------------------------------------

Double_t omega_cm(Double_t omega_lab)
{
  Double_t s = MASS2_INITIAL + 2.0*omega_lab*MASS_INITIAL;
  Double_t sqrt_s = Sqrt(s);
  Double_t omega = (s - MASS2_INITIAL)/(2.0*sqrt_s);

  return omega;
}

//-----------------------------------------------------------------------------

Double_t omega_lab(Double_t W_cm)
{
  Double_t s = W_cm*W_cm;
  Double_t omega = (s - MASS2_INITIAL)/(2.0*MASS_INITIAL);

  return omega;
}

//-----------------------------------------------------------------------------

Double_t rho(Double_t Omega)
{
  Double_t q = q_mes(Omega);
  Double_t k = omega_cm(Omega);

  return (q/k);
}

//-----------------------------------------------------------------------------

Double_t L(Int_t l, Double_t x)
{
  if(l==0) return 1.0; //Legendre polynom P_0
  if(l==1) return x;   //Legendre polynom P_1
  if(l > 1) return (1.0/l) *((2.0*l - 1.0)*x*L(l-1, x) - (l - 1.0)*L(l-2, x));
  return 0.0; //Any other case, e.g. for l < 0
}

//-----------------------------------------------------------------------------

Double_t DL(Int_t l, Double_t x)
{
  if(l==0) return 0.0; //Legendre polynom 1st derivative P'_0
  if(l==1) return 1.0; //Legendre polynom 1st derivative P'_1
  if(l > 1) return (1.0/(l - 1.0)) * ((2.0*l - 1.0)*x*DL(l-1, x) - l*DL(l-2, x));
  return 0.0; //Any other case, e.g. for l < 0
}

//-----------------------------------------------------------------------------

Double_t D2L(Int_t l, Double_t x)
{
  if(l==0) return 0.0; //Legendre polynom 2nd derivative P"_0
  if(l==1) return 0.0; //Legendre polynom 2nd derivative P"_1
  if(l==2) return 3.0; //Legendre polynom 2nd derivative P"_2
  if(l > 2) return (1.0/(l - 2.0)) * ((2.0*l - 1.0)*x*D2L(l-1, x) - (l + 1.0)*D2L(l-2, x));
  return 0.0; //Any other case, e.g. for l < 0
}

//-----------------------------------------------------------------------------

//MAID CGLN amplitude F1 in expansion up to l_max
TComplex F1(Double_t CosTheta, Int_t e)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=0; l<L_MAX+1; l++)
    Complex+=((maid_Mp[l][e]*(1.0*l) + maid_Ep[l][e]) * DL(l+1, CosTheta) + (maid_Mm[l][e]*(1.0*l + 1.0) + maid_Em[l][e]) * DL(l-1, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//MAID CGLN amplitude F2 in expansion up to l_max
TComplex F2(Double_t CosTheta, Int_t e)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=1; l<L_MAX+1; l++)
    Complex+=((maid_Mp[l][e]*(1.0*l + 1.0) + maid_Mm[l][e]*(1.0*l)) * DL(l, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//MAID CGLN amplitude F3 in expansion up to l_max
TComplex F3(Double_t CosTheta, Int_t e)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=1; l<L_MAX+1; l++)
    Complex+=((maid_Ep[l][e] - maid_Mp[l][e]) * D2L(l+1, CosTheta) + (maid_Em[l][e] + maid_Mm[l][e]) * D2L(l-1, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

//MAID CGLN amplitude F4 in expansion up to l_max
TComplex F4(Double_t CosTheta, Int_t e)
{
  TComplex Complex(0.0, 0.0);

  for(Int_t l=2; l<L_MAX+1; l++)
    Complex+=((maid_Mp[l][e] - maid_Ep[l][e] - maid_Mm[l][e] - maid_Em[l][e]) * D2L(l, CosTheta));

  return Complex;
}

//-----------------------------------------------------------------------------

TComplex F1cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(F1(CosTheta, e)); }
TComplex F2cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(F2(CosTheta, e)); }
TComplex F3cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(F3(CosTheta, e)); }
TComplex F4cc(Double_t CosTheta, Int_t e){ return TComplex::Conjugate(F4(CosTheta, e)); }

//-----------------------------------------------------------------------------

void Parse_MAID()
{
  Char_t Buffer[1024];
  Double_t Energy, Re, Im;
  FILE* MAID_Ep[LBINS];
  FILE* MAID_Em[LBINS];
  FILE* MAID_Mp[LBINS];
  FILE* MAID_Mm[LBINS];

  printf("Loading model multipoles... ");

  //Open s,p waves
  MAID_Ep[0] = fopen("model/E0p.txt", "r");
  MAID_Ep[1] = fopen("model/E1p.txt", "r");
  MAID_Mp[1] = fopen("model/M1p.txt", "r");
  MAID_Mm[1] = fopen("model/M1m.txt", "r");
  //Open d,f,g,... waves
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    sprintf(Buffer, "model/E%dp.txt", l); MAID_Ep[l] = fopen(Buffer, "r");
    sprintf(Buffer, "model/E%dm.txt", l); MAID_Em[l] = fopen(Buffer, "r");
    sprintf(Buffer, "model/M%dp.txt", l); MAID_Mp[l] = fopen(Buffer, "r");
    sprintf(Buffer, "model/M%dm.txt", l); MAID_Mm[l] = fopen(Buffer, "r");
  }

  //Load E0+ multipole
  maid_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), MAID_Ep[0]);
  fgets(Buffer, sizeof(Buffer), MAID_Ep[0]);
  while(!feof(MAID_Ep[0]))
  {
    if(fscanf(MAID_Ep[0], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
    maid_Ep[0][maid_bin] = TComplex(Re, Im);
    maid_en[maid_bin] = omega_lab(Energy);
    maid_bin++;
  }

  //Load E1+ multipole
  maid_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), MAID_Ep[1]);
  fgets(Buffer, sizeof(Buffer), MAID_Ep[1]);
  while(!feof(MAID_Ep[1]))
  {
    if(fscanf(MAID_Ep[1], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
    maid_Ep[1][maid_bin] = TComplex(Re, Im);
    maid_en[maid_bin] = omega_lab(Energy);
    maid_bin++;
  }

  //Load M1+ multipole
  maid_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), MAID_Mp[1]);
  fgets(Buffer, sizeof(Buffer), MAID_Mp[1]);
  while(!feof(MAID_Mp[1]))
  {
    if(fscanf(MAID_Mp[1], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
    maid_Mp[1][maid_bin] = TComplex(Re, Im);
    maid_en[maid_bin] = omega_lab(Energy);
    maid_bin++;
  }

  //Load M1- multipole
  maid_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), MAID_Mm[1]);
  fgets(Buffer, sizeof(Buffer), MAID_Mm[1]);
  while(!feof(MAID_Mm[1]))
  {
    if(fscanf(MAID_Mm[1], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
    maid_Mm[1][maid_bin] = TComplex(Re, Im);
    maid_en[maid_bin] = omega_lab(Energy);
    maid_bin++;
  }

   //Load El+ multipole for d,f,... waves
  for(Int_t l=2;l<L_MAX+1; l++)
  {
    maid_bin = 0;
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), MAID_Ep[l]);
    fgets(Buffer, sizeof(Buffer), MAID_Ep[l]);
    while(!feof(MAID_Ep[l]))
    {
      if(fscanf(MAID_Ep[l], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
      maid_Ep[l][maid_bin] = TComplex(Re, Im);
      maid_en[maid_bin] = omega_lab(Energy);
      maid_bin++;
    }
  }

   //Load El- multipole for d,f,... waves
  for(Int_t l=2;l<L_MAX+1; l++)
  {
    maid_bin = 0;
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), MAID_Em[l]);
    fgets(Buffer, sizeof(Buffer), MAID_Em[l]);
    while(!feof(MAID_Em[l]))
    {
      if(fscanf(MAID_Em[l], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
      maid_Em[l][maid_bin] = TComplex(Re, Im);
      maid_en[maid_bin] = omega_lab(Energy);
      maid_bin++;
    }
  }

   //Load Ml+ multipole for d,f,... waves
  for(Int_t l=2;l<L_MAX+1; l++)
  {
    maid_bin = 0;
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), MAID_Mp[l]);
    fgets(Buffer, sizeof(Buffer), MAID_Mp[l]);
    while(!feof(MAID_Mp[l]))
    {
      if(fscanf(MAID_Mp[l], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
      maid_Mp[l][maid_bin] = TComplex(Re, Im);
      maid_en[maid_bin] = omega_lab(Energy);
      maid_bin++;
    }
  }

   //Load Ml- multipole for d,f,... waves
  for(Int_t l=2;l<L_MAX+1; l++)
  {
    maid_bin = 0;
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), MAID_Mm[l]);
    fgets(Buffer, sizeof(Buffer), MAID_Mm[l]);
    while(!feof(MAID_Mm[l]))
    {
      if(fscanf(MAID_Mm[l], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
      maid_Mm[l][maid_bin] = TComplex(Re, Im);
      maid_en[maid_bin] = omega_lab(Energy);
      maid_bin++;
    }
  }

  //Close s,p waves
  fclose(MAID_Ep[0]);
  fclose(MAID_Ep[1]);
  fclose(MAID_Mp[1]);
  fclose(MAID_Mm[1]);
  //Close d,f,... waves
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    fclose(MAID_Ep[l]);
    fclose(MAID_Em[l]);
    fclose(MAID_Mp[l]);
    fclose(MAID_Mm[l]);
  }

  printf("%2d multipoles at %4d energies loaded\n", 4*L_MAX, maid_bin);
  return;

  //Debug output
  printf("EBins: %d\n", maid_bin);
  for(Int_t e=0; e<maid_bin; e++)
  {
    printf("%d (%f MeV)\n", e, maid_en[e]);
    printf("%f %f %f %f %f %f %f %f\n",
           maid_Ep[0][e].Re(), maid_Ep[1][e].Re(), maid_Mp[1][e].Re(), maid_Mm[1][e].Re(), maid_Ep[2][e].Re(), maid_Mp[2][e].Re(), maid_Em[2][e].Re(), maid_Mm[2][e].Re());
    printf("%f %f %f %f %f %f %f %f\n",
           maid_Ep[0][e].Im(), maid_Ep[1][e].Im(), maid_Mp[1][e].Im(), maid_Mm[1][e].Im(), maid_Ep[2][e].Im(), maid_Mp[2][e].Im(), maid_Em[2][e].Im(), maid_Mm[2][e].Im());
  }
}

//-----------------------------------------------------------------------------

Double_t sigma0(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F1(CosTheta, e).Rho2() + F2(CosTheta, e).Rho2()
                   + Sin2Theta*(0.5*F3(CosTheta, e).Rho2() + 0.5*F4(CosTheta, e).Rho2() + F2cc(CosTheta, e)*F3(CosTheta, e) + F1cc(CosTheta, e)*F4(CosTheta, e) + CosTheta*F3cc(CosTheta, e)*F4(CosTheta, e))
                   - 2.0*CosTheta*F1cc(CosTheta, e)*F2(CosTheta, e);

  return Complex.Re()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------



Double_t sigmaS(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = 0.5*(F3(CosTheta, e).Rho2() + F4(CosTheta, e).Rho2()) + F2cc(CosTheta, e)*F3(CosTheta, e) + F1cc(CosTheta, e)*F4(CosTheta, e) + CosTheta*F3cc(CosTheta, e)*F4(CosTheta, e);
  return -Sin2Theta*Complex.Re()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaT(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F1cc(CosTheta, e)*F3(CosTheta, e) - F2cc(CosTheta, e)*F4(CosTheta, e)
                   + CosTheta*(F1cc(CosTheta, e)*F4(CosTheta, e) - F2cc(CosTheta, e)*F3(CosTheta, e))
                   - Sin2Theta*F3cc(CosTheta, e)*F4(CosTheta, e);
  return SinTheta*Complex.Im()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaP(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = 2.0*F1cc(CosTheta, e)*F2(CosTheta, e) + F1cc(CosTheta, e)*F3(CosTheta, e) - F2cc(CosTheta, e)*F4(CosTheta, e)
                   - CosTheta*(F2cc(CosTheta, e)*F3(CosTheta, e) - F1cc(CosTheta, e)*F4(CosTheta, e))
                   - Sin2Theta*F3cc(CosTheta, e)*F4(CosTheta, e);
  return -SinTheta*Complex.Im()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaH(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());

  TComplex Complex = 2.0*F1cc(CosTheta, e)*F2(CosTheta, e) + F1cc(CosTheta, e)*F3(CosTheta, e) - F2cc(CosTheta, e)*F4(CosTheta, e)
                   + CosTheta*(F1cc(CosTheta, e)*F4(CosTheta, e) - F2cc(CosTheta, e)*F3(CosTheta, e));
  return SinTheta*Complex.Im()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaG(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F2cc(CosTheta, e)*F3(CosTheta, e) + F1cc(CosTheta, e)*F4(CosTheta, e);
  return Sin2Theta*Complex.Im()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaF(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());

  TComplex Complex = F1cc(CosTheta, e)*F3(CosTheta, e) - F2cc(CosTheta, e)*F4(CosTheta, e)
                   - CosTheta*(F2cc(CosTheta, e)*F3(CosTheta, e) - F1cc(CosTheta, e)*F4(CosTheta, e));
  return SinTheta*Complex.Re()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaE(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F1(CosTheta, e).Rho2() + F2(CosTheta, e).Rho2()
                   - 2.0*CosTheta*F1cc(CosTheta, e)*F2(CosTheta, e)
                   + Sin2Theta*(F2cc(CosTheta, e)*F3(CosTheta, e) + F1cc(CosTheta, e)*F4(CosTheta, e));
  return Complex.Re()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaCx(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());

  TComplex Complex = F1(CosTheta, e).Rho2() - F2(CosTheta, e).Rho2()
                   + F1cc(CosTheta, e)*F4(CosTheta, e) - F2cc(CosTheta, e)*F3(CosTheta, e)
                   + CosTheta*(F1cc(CosTheta, e)*F3(CosTheta, e) - F2cc(CosTheta, e)*F4(CosTheta, e));
  return SinTheta*Complex.Re()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaCz(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = 2.0*F1cc(CosTheta, e)*F2(CosTheta, e)
                   + Sin2Theta*(F1cc(CosTheta, e)*F3(CosTheta, e) + F2cc(CosTheta, e)*F4(CosTheta, e))
                   - CosTheta*(F1(CosTheta, e).Rho2() + F2(CosTheta, e).Rho2());
  return Complex.Re()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaOx(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());

  TComplex Complex = F1cc(CosTheta, e)*F4(CosTheta, e) - F2cc(CosTheta, e)*F3(CosTheta, e)
                   + CosTheta*(F1cc(CosTheta, e)*F3(CosTheta, e) - F2cc(CosTheta, e)*F4(CosTheta, e));
  return -SinTheta*Complex.Im()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t sigmaOz(Double_t Theta, Int_t e)
{
  Double_t SinTheta = Sin(Theta*DegToRad());
  Double_t CosTheta = Cos(Theta*DegToRad());
  Double_t Sin2Theta = SinTheta*SinTheta;

  TComplex Complex = F1cc(CosTheta, e)*F3(CosTheta, e) + F2cc(CosTheta, e)*F4(CosTheta, e);
  return -Sin2Theta*Complex.Im()*rho(maid_en[e])*UNIT;
}

//-----------------------------------------------------------------------------

Double_t S(Double_t Theta, Int_t e)
{
  return sigmaS(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t T(Double_t Theta, Int_t e)
{
  return sigmaT(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t P(Double_t Theta, Int_t e)
{
  return sigmaP(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t H(Double_t Theta, Int_t e)
{
  return sigmaH(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t G(Double_t Theta, Int_t e)
{
  return sigmaG(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t F(Double_t Theta, Int_t e)
{
  return sigmaF(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t E(Double_t Theta, Int_t e)
{
  return sigmaE(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t Cx(Double_t Theta, Int_t e)
{
  return sigmaCx(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t Cz(Double_t Theta, Int_t e)
{
  return sigmaCz(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t Ox(Double_t Theta, Int_t e)
{
  return sigmaOx(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

Double_t Oz(Double_t Theta, Int_t e)
{
  return sigmaOz(Theta, e)/sigma0(Theta, e);
}

//-----------------------------------------------------------------------------

void Print_sg0()
{
  FILE* Out = fopen("sg0.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = sigma0(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_S()
{
  FILE* Out = fopen("S.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = S(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_T()
{
  FILE* Out = fopen("T.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = T(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_P()
{
  FILE* Out = fopen("P.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = P(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_E()
{
  FILE* Out = fopen("E.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = E(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_F()
{
  FILE* Out = fopen("F.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = F(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_G()
{
  FILE* Out = fopen("G.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = G(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_H()
{
  FILE* Out = fopen("H.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = H(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_Cx()
{
  FILE* Out = fopen("Cx.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = Cx(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_Cz()
{
  FILE* Out = fopen("Cz.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = Cz(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_Ox()
{
  FILE* Out = fopen("Ox.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = Ox(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------

void Print_Oz()
{
  FILE* Out = fopen("Oz.txt", "w");
  Double_t th;
  Double_t Obs, Obs_meas;

  Parse_MAID();

  for(Int_t e=0; e<maid_bin; e++)
  {
    fprintf(Out, "E = %8.3f MeV, Wght = 1.00, Syst = 0.05\n", maid_en[e]);
    for(Int_t t=0; t<18; t++)
    {
      th = 10.0*t + 5.0;
      Obs  = Oz(th, e);
      Obs_meas = Obs*(1.0 + gRandom->Gaus(0.0, 0.05/Sin(th*DegToRad())));
      fprintf(Out, "%7.3f  %10.6f  %10.6f\n", th, Obs_meas, fabs(Obs)*0.05/Sin(th*DegToRad()));
    }
    fprintf(Out, "------------------------------------------\n");
  }
  fclose(Out);
}

//-----------------------------------------------------------------------------
