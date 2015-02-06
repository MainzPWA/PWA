static const Int_t N_MAX = 1024;

//-----------------------------------------------------------------------------

Double_t ReA_ppi0[N_MAX];
Double_t ImA_ppi0[N_MAX];
Double_t ReA_npip[N_MAX];
Double_t ImA_npip[N_MAX];
Double_t ReDA_ppi0[N_MAX];
Double_t ImDA_ppi0[N_MAX];
Double_t ReDA_npip[N_MAX];
Double_t ImDA_npip[N_MAX];
Double_t W_ppi0[N_MAX];
Double_t W_npip[N_MAX];
Int_t Pts_npip;
Int_t Pts_ppi0;

Double_t ReA_12[N_MAX];
Double_t ImA_12[N_MAX];
Double_t ReA_32[N_MAX];
Double_t ImA_32[N_MAX];
Double_t ReDA_12[N_MAX];
Double_t ImDA_12[N_MAX];
Double_t ReDA_32[N_MAX];
Double_t ImDA_32[N_MAX];
Double_t W_iso[N_MAX];
Int_t Pts_iso;

//-----------------------------------------------------------------------------

void IsoExtract(Int_t L_MAX = 5, Double_t W_LO=0.0, Double_t W_HI=2500.0)
{
  Char_t Buffer[256];

  Isospin("E0p", W_LO, W_HI);
  Isospin("E1p", W_LO, W_HI);
  Isospin("M1m", W_LO, W_HI);
  Isospin("M1p", W_LO, W_HI);

  for(Int_t l=2; l<L_MAX+1; l++)
  {
    sprintf(Buffer, "E%1dm", l);
    Isospin(Buffer, W_LO, W_HI);
    sprintf(Buffer, "E%1dp", l);
    Isospin(Buffer, W_LO, W_HI);
    sprintf(Buffer, "M%1dm", l);
    Isospin(Buffer, W_LO, W_HI);
    sprintf(Buffer, "M%1dp", l);
    Isospin(Buffer, W_LO, W_HI);
  }
}

//-----------------------------------------------------------------------------

void Isospin(Char_t* Mlp, Double_t W_LO, Double_t W_HI)
{
  Char_t Buffer[256];
  Char_t Fancy[16];
  Int_t n_0, n_p;
  FILE* Out_12;
  FILE* Out_32;

  //Open multipole files for p pi0 and n pi+ channels
  Load_ppi0(Mlp);
  Load_npip(Mlp);

  Pts_iso = 0;
  //If p pi0 channel contains more data (i.e. has finer energy binning)...
  if(Pts_ppi0 > Pts_npip)
  {
    //...perform multipole extraction at p pi0 energies
    for(n_0=0; n_0<Pts_ppi0; n_0++)
    {
      //Check requested energy bounds
      if((W_ppi0[n_0] < W_LO) || (W_ppi0[n_0] > W_HI)) continue;
      //Find n pi+ entry with closest energy
      n_p = GetBin_npip(W_ppi0[n_0]);
      //Get isospin multipoles...
      ReA_12[Pts_iso] = A_12(ReA_ppi0[n_0], ReA_npip[n_p]);
      ReA_32[Pts_iso] = A_32(ReA_ppi0[n_0], ReA_npip[n_p]);
      ImA_12[Pts_iso] = A_12(ImA_ppi0[n_0], ImA_npip[n_p]);
      ImA_32[Pts_iso] = A_32(ImA_ppi0[n_0], ImA_npip[n_p]);
      //...and multipole errors
      ReDA_12[Pts_iso] = DA_12(ReDA_ppi0[n_0], ReDA_npip[n_p]);
      ReDA_32[Pts_iso] = DA_32(ReDA_ppi0[n_0], ReDA_npip[n_p]);
      ImDA_12[Pts_iso] = DA_12(ImDA_ppi0[n_0], ImDA_npip[n_p]);
      ImDA_32[Pts_iso] = DA_32(ImDA_ppi0[n_0], ImDA_npip[n_p]);
      //Store energy and count datapoints
      W_iso[Pts_iso] = W_ppi0[n_0];
      Pts_iso++;
    }
  }
  //If n pi+ channel contains more data (i.e. has finer energy binning)...
  else if(Pts_ppi0 < Pts_npip)
  {
    //...perform multipole extraction at n pi+ energies
    for(n_p=0; n_p<Pts_npip; n_p++)
    {
      //Check requested energy bounds
      if((W_npip[n_p] < W_LO) || (W_npip[n_p] > W_HI)) continue;
      //Find p pi0 entry with closest energy
      n_0 = GetBin_ppi0(W_npip[n_p]);
      //Get isospin multipoles...
      ReA_12[Pts_iso] = A_12(ReA_ppi0[n_0], ReA_npip[n_p]);
      ReA_32[Pts_iso] = A_32(ReA_ppi0[n_0], ReA_npip[n_p]);
      ImA_12[Pts_iso] = A_12(ImA_ppi0[n_0], ImA_npip[n_p]);
      ImA_32[Pts_iso] = A_32(ImA_ppi0[n_0], ImA_npip[n_p]);
      //...and multipole errors
      ReDA_12[Pts_iso] = DA_12(ReDA_ppi0[n_0], ReDA_npip[n_p]);
      ReDA_32[Pts_iso] = DA_32(ReDA_ppi0[n_0], ReDA_npip[n_p]);
      ImDA_12[Pts_iso] = DA_12(ImDA_ppi0[n_0], ImDA_npip[n_p]);
      ImDA_32[Pts_iso] = DA_32(ImDA_ppi0[n_0], ImDA_npip[n_p]);
      //Store energy and count datapoints
      W_iso[Pts_iso] = W_npip[n_p];
      Pts_iso++;
    }
  }
  //Sort(0, Pts_iso-1);

  //Create nicer multipole name (e.g. 'E0+' from 'E0p')
  strcpy(Fancy, Mlp);
  for(Int_t t=0; t<strlen(Fancy); t++)
  {
    if(Fancy[t]=='p') Fancy[t] = '+';
    if(Fancy[t]=='m') Fancy[t] = '-';
  }

  //Open output files for p1/2 and 3/2 multipoles
  sprintf(Buffer, "isospin/%s_p12.txt", Mlp);
  Out_12 = fopen(Buffer, "w");
  sprintf(Buffer, "isospin/%s_32.txt", Mlp);
  Out_32 = fopen(Buffer, "w");

  //Print table headers for p1/2 and 3/2 multipoles
  fprintf(Out_12, "   W                 %s_p1/2\n", Fancy);
  fprintf(Out_12, " (MeV)      Re       DRe      Im       DIm\n");
  fprintf(Out_32, "   W                 %s_3/2\n", Fancy);
  fprintf(Out_32, " (MeV)      Re       DRe      Im       DIm\n");

  //Print p1/2 and 3/2 multipole values and errors to output files
  for(Int_t n=0; n<Pts_iso; n++)
  {
    fprintf(Out_12, "%8.3f  %7.3f  %7.3f  %7.3f  %7.3f\n", W_iso[n], ReA_12[n], ReDA_12[n], ImA_12[n], ImDA_12[n]);
    fprintf(Out_32, "%8.3f  %7.3f  %7.3f  %7.3f  %7.3f\n", W_iso[n], ReA_32[n], ReDA_32[n], ImA_32[n], ImDA_32[n]);
  }

  //Close output files for p1/2 and 3/2 multipoles
  fclose(Out_12);
  fclose(Out_32);
}

//-----------------------------------------------------------------------------

void IsoMultipole(Char_t* Mlp, Char_t* Iso, Bool_t SAVE=false, Double_t Lo=0.0, Double_t Hi=0.0)
{
  if(SAVE) TCanvas* Canvas = new TCanvas();
  //if(SAVE) SetBit(TH1::kNoTitle);

  FILE* InPlots;
  FILE* InModel_0;
  FILE* InModel_p;
  Char_t Buffer[256];
  Char_t M, p;
  Int_t l;
  Double_t PlRe[N_MAX];
  Double_t PlIm[N_MAX];
  Double_t PlDRe[N_MAX];
  Double_t PlDIm[N_MAX];
  Double_t PlW[N_MAX];
  Double_t MoRe_0[N_MAX];
  Double_t MoRe_p[N_MAX];
  Double_t MoIm_0[N_MAX];
  Double_t MoIm_p[N_MAX];
  Double_t MoW_0[N_MAX];
  Double_t MoW_p[N_MAX];
  Double_t MoRe[N_MAX];
  Double_t MoIm[N_MAX];
  Int_t PlPts;
  Int_t MoPts_0, MoPts_p;
  Double_t W, Re, DRe, Im, DIm;
  Double_t Min = 0.0;
  Double_t Max = 0.0;
  TGraphErrors* PlotsRe;
  TGraphErrors* PlotsIm;
  TGraph* ModelRe;
  TGraph* ModelIm;

  //Decompose multipole name
  sscanf(Mlp, "%c%d%c", &M, &l, &p);

  //Open text file with fit results for given multipole
  sprintf(Buffer, "isospin/%s_%s.txt", Mlp, Iso);
  InPlots = fopen(Buffer, "r");

  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), InPlots);
  fgets(Buffer, sizeof(Buffer), InPlots);
  //Read multipole fit values from file
  PlPts = 0;
  while(!feof(InPlots))
  {
    if(fscanf(InPlots, "%lf %lf %lf %lf %lf", &W, &Re, &DRe, &Im, &DIm)==5)
    {
      //Look for maximum & minimum of values to adjust plotting range
      if(Re+DRe > Max) Max = Re+DRe;
      if(Im+DIm > Max) Max = Im+DIm;
      if(Re-DRe < Min) Min = Re-DRe;
      if(Im-DIm < Min) Min = Im-DIm;
      //Add real and imaginary parts of multipole (w/ errors) to graph
      PlW[PlPts] = W;
      PlRe[PlPts] = Re;
      PlIm[PlPts] = Im;
      PlDRe[PlPts] = DRe;
      PlDIm[PlPts] = DIm;
      PlPts++;
    }
  }
  //Close file with fit values
  fclose(InPlots);

  //Create graphs for real and imaginary parts of fitted multipole
  PlotsRe = new TGraphErrors(PlPts, PlW, PlRe, NULL, PlDRe);
  PlotsIm = new TGraphErrors(PlPts, PlW, PlIm, NULL, PlDIm);

  //Color, line size, marker style adjustments
  PlotsRe->SetLineColor(kRed+2);
  PlotsIm->SetLineColor(kBlue+2);
  PlotsRe->SetMarkerColor(kRed+2);
  PlotsIm->SetMarkerColor(kBlue+2);
  PlotsRe->SetMarkerStyle(21);
  PlotsIm->SetMarkerStyle(21);
  PlotsRe->SetMarkerSize(0.7);
  PlotsIm->SetMarkerSize(0.7);

  //Set plot titles and object names
  sprintf(Buffer, "%s_%s", Mlp, Iso);
  PlotsRe->SetTitle(Buffer);
  PlotsIm->SetTitle(Buffer);
  sprintf(Buffer, "Fit_Re%s_%s", Mlp, Iso);
  PlotsRe->SetName(Buffer);
  sprintf(Buffer, "Fit_Im%s_%s", Mlp, Iso);
  PlotsIm->SetName(Buffer);

  PlotsRe->Draw("APZ"); //Plot with x,y axis, points and small error bars
  PlotsIm->Draw("PZ"); //Plot with points and small error bars, into existing frame

  //x-axis labeling
  PlotsRe->GetXaxis()->SetTitle("W / MeV");
  //y-axis labeling
  sprintf(Buffer, "%c_{%1d%c}", M, l, p);
  for(Int_t n=0; n<strlen(Buffer); n++)
  {
    if(Buffer[n]=='p') Buffer[n] = '+';
    if(Buffer[n]=='m') Buffer[n] = '-';
  }
  if(!strcmp(Iso, "32"))  sprintf(Buffer, "%s^{3/2}",  Buffer, Iso);
  if(!strcmp(Iso, "p12")) sprintf(Buffer, "%s^{p1/2}", Buffer, Iso);
  PlotsRe->GetYaxis()->SetTitle(Buffer);

  //Open text file with model values for given p pi0 multipole
  sprintf(Buffer, "model/ppi0/%s.txt", Mlp);
  InModel_0 = fopen(Buffer, "r");
  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), InModel_0);
  fgets(Buffer, sizeof(Buffer), InModel_0);
  //Read multipole model values from file
  MoPts_0 = 0;
  while(!feof(InModel_0))
  {
    if(fscanf(InModel_0, "%lf %lf %lf", &W, &Re, &Im)==3)
    {
      MoW_0[MoPts_0]  = W;
      MoRe_0[MoPts_0] = Re;
      MoIm_0[MoPts_0] = Im;
      MoPts_0++;
    }
  }
  //Close file with model values
  fclose(InModel_0);

  //Open text file with model values for given n pi+ multipole
  sprintf(Buffer, "model/npip/%s.txt", Mlp);
  InModel_p = fopen(Buffer, "r");
  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), InModel_p);
  fgets(Buffer, sizeof(Buffer), InModel_p);
  //Read multipole model values from file
  MoPts_p = 0;
  while(!feof(InModel_p))
  {
    if(fscanf(InModel_p, "%lf %lf %lf", &W, &Re, &Im)==3)
    {
      MoW_p[MoPts_p]  = W;
      MoRe_p[MoPts_p] = Re;
      MoIm_p[MoPts_p] = Im;
      MoPts_p++;
    }
  }
  //Close file with model values
  fclose(InModel_p);

  //Create selected isospin multipole from p pi0 and n pi+ multipoles
  for(Int_t wp=0; wp<MoPts_p; wp++)
  {
    //Find corresponding energy bin between n pi+ and p pi0 multipoles
    Int_t w0;
    for(Int_t w0=0; w0<MoPts_0; w0++)
      if(MoW_p[wp]==MoW_0[w0]) break;
    //Create isospin multipoles
    if(!strcmp(Iso, "32"))  { MoRe[wp] = A_32(MoRe_0[w0], MoRe_p[wp]); MoIm[wp] = A_32(MoIm_0[w0], MoIm_p[wp]); }
    if(!strcmp(Iso, "p12")) { MoRe[wp] = A_12(MoRe_0[w0], MoRe_p[wp]); MoIm[wp] = A_12(MoIm_0[w0], MoIm_p[wp]); }
  }

  //Create graphs for real and imaginary parts of model multipole
  ModelRe = new TGraph(MoPts_p, MoW_p, MoRe);
  ModelIm = new TGraph(MoPts_p, MoW_p, MoIm);

  //Color, line size adjustments
  ModelRe->SetLineColor(kRed);
  ModelIm->SetLineColor(kBlue);
  ModelRe->SetLineWidth(2);
  ModelIm->SetLineWidth(2);

  //Set plot titles and object names
  sprintf(Buffer, "%s_%s", Mlp, Iso);
  ModelRe->SetTitle(Buffer);
  ModelIm->SetTitle(Buffer);
  sprintf(Buffer, "Model_Re%s_%s", Mlp, Iso);
  ModelRe->SetName(Buffer);
  sprintf(Buffer, "Model_Im%s_%s", Mlp, Iso);
  ModelIm->SetName(Buffer);

  //Plot graphs
  ModelRe->Draw("L"); //Plot as line, into same frame
  ModelIm->Draw("L"); //Plot as line, into same frame

  //Adjust drawing ranges for y-axis
  if((Lo==0.0) && (Hi==0.0))
  {
    PlotsRe->SetMinimum(Min*1.05);
    PlotsRe->SetMaximum(Max*1.05);
  }
  else
  {
    PlotsRe->SetMinimum(Lo);
    PlotsRe->SetMaximum(Hi);
  }

  if(SAVE) PlotsRe->SetTitle("");
  if(SAVE) sprintf(Buffer, "isospin/%s_%s.eps", Mlp, Iso);
  if(SAVE) Canvas->SaveAs(Buffer);
}

//-----------------------------------------------------------------------------

void Load_ppi0(Char_t* Mlp)
{
  Char_t Buffer[256];
  Double_t W, Re, DRe, Im, DIm;
  FILE* MlpFile;

  Pts_ppi0 = 0;
  sprintf(Buffer, "data/ppi0/%s.txt", Mlp);
  MlpFile = fopen(Buffer, "r");

  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), MlpFile);
  fgets(Buffer, sizeof(Buffer), MlpFile);
  //Read multipole fit values from file
  while(!feof(MlpFile))
  {
    if(fscanf(MlpFile, "%lf %lf %lf %lf %lf", &W, &Re, &DRe, &Im, &DIm)==5)
    {
      W_ppi0[Pts_ppi0] = W;
      ReA_ppi0[Pts_ppi0] = Re;
      ImA_ppi0[Pts_ppi0] = Im;
      ReDA_ppi0[Pts_ppi0] = DRe;
      ImDA_ppi0[Pts_ppi0] = DIm;
      Pts_ppi0++;
    }
  }
  fclose(MlpFile);
}

//-----------------------------------------------------------------------------

void Load_npip(Char_t* Mlp)
{
  Char_t Buffer[256];
  Double_t W, Re, DRe, Im, DIm;
  FILE* MlpFile;

  Pts_npip = 0;
  sprintf(Buffer, "data/npip/%s.txt", Mlp);
  MlpFile = fopen(Buffer, "r");

  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), MlpFile);
  fgets(Buffer, sizeof(Buffer), MlpFile);
  //Read multipole fit values from file
  while(!feof(MlpFile))
  {
    if(fscanf(MlpFile, "%lf %lf %lf %lf %lf", &W, &Re, &DRe, &Im, &DIm)==5)
    {
      W_npip[Pts_npip] = W;
      ReA_npip[Pts_npip] = Re;
      ImA_npip[Pts_npip] = Im;
      ReDA_npip[Pts_npip] = DRe;
      ImDA_npip[Pts_npip] = DIm;
      Pts_npip++;
    }
  }
  fclose(MlpFile);
}

//-----------------------------------------------------------------------------

Int_t GetBin_npip(Double_t W)
{
  Double_t Min = 1e38;
  Int_t Bin = 0;

  for(Int_t e=0; e<Pts_npip; e++)
    if(fabs(W_npip[e] - W) < Min)
    {
      Min = fabs(W_npip[e] - W);
      Bin = e;
    }
  return Bin;
}

//-----------------------------------------------------------------------------

Int_t GetBin_ppi0(Double_t W)
{
  Double_t Min = 1e38;
  Int_t Bin = 0;

  for(Int_t e=0; e<Pts_ppi0; e++)
    if(fabs(W_ppi0[e] - W) < Min)
    {
      Min = fabs(W_ppi0[e] - W);
      Bin = e;
    }
  return Bin;
}

//-----------------------------------------------------------------------------

Double_t A_12(Double_t A_ppi0, Double_t A_npip)
{
  return (1.0/3.0)*(TMath::Sqrt(2.0)*A_npip + A_ppi0);
}

//-----------------------------------------------------------------------------

Double_t DA_12(Double_t DA_ppi0, Double_t DA_npip)
{
  Double_t npip = TMath::Sqrt(2.0)*DA_npip;
  Double_t ppi0 = DA_ppi0;

  return (1.0/3.0)*TMath::Sqrt(npip*npip + ppi0*ppi0);
}

//-----------------------------------------------------------------------------

Double_t A_32(Double_t A_ppi0, Double_t A_npip)
{
  return A_ppi0 - TMath::Sqrt(0.5)*A_npip;
}

//-----------------------------------------------------------------------------

Double_t DA_32(Double_t DA_ppi0, Double_t DA_npip)
{
  Double_t ppi0 = DA_ppi0;
  Double_t npip = TMath::Sqrt(0.5)*DA_npip;

  return TMath::Sqrt(ppi0*ppi0 + npip*npip);
}

//-----------------------------------------------------------------------------

//void Sort(Int_t l, Int_t r) //Quicksort implementation
//{
//  if(r > l)
//  {
//    Int_t i = l-1;
//    Int_t j = r;
//
//   for(;;)
//   {
//     while(W_iso[++i] < W_iso[r]);
//     while((W_iso[--j] > W_iso[r]) && (j>i));
//     if(i>=j) break;
//     Swap(&W_iso[i], &W_iso[j]);
//     Swap(&ReA_12[i], &ReA_12[j]);
//     Swap(&ImA_12[i], &ImA_12[j]);
//     Swap(&ReA_32[i], &ReA_32[j]);
//     Swap(&ImA_32[i], &ImA_32[j]);
//     Swap(&ReDA_12[i], &ReDA_12[j]);
//     Swap(&ImDA_12[i], &ImDA_12[j]);
//     Swap(&ReDA_32[i], &ReDA_32[j]);
//     Swap(&ImDA_32[i], &ImDA_32[j]);
//   }
//   Swap(&W_iso[i],  &W_iso[r]);
//   Swap(&ReA_12[i], &ReA_12[r]);
//   Swap(&ImA_12[i], &ImA_12[r]);
//   Swap(&ReA_32[i], &ReA_32[r]);
//   Swap(&ImA_32[i], &ImA_32[r]);
//   Swap(&ReDA_12[i], &ReDA_12[r]);
//   Swap(&ImDA_12[i], &ImDA_12[r]);
//   Swap(&ReDA_32[i], &ReDA_32[r]);
//   Swap(&ImDA_32[i], &ImDA_32[r]);
//
//   Sort(l, i-1);
//   Sort(i+1, r);
//  }
//}

//-----------------------------------------------------------------------------

//void Swap(Double_t* a, Double_t* b)
//{
//  Double_t tmp;
//
//  tmp = *a;
//  *a = *b;
//  *b = tmp;
//}

//-----------------------------------------------------------------------------
