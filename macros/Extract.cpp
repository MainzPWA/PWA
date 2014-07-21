static const Double_t MASS_PROTON  = 938.2720;
static const Double_t MASS_NEUTRON = 939.5654;
static const Int_t N_MAX = 1024;
static const Int_t LBINS = 10;
static const Int_t SOL = 16;

static const Int_t RED[SOL]  = { kRed+2,  kRed-3,  kRed-6,  kRed-9,  kRed-10, kRed-10, kRed-10, kRed-10,
                                 kRed-10, kRed-10, kRed-10, kRed-10, kRed-10, kRed-10, kRed-10, kRed-10 };
static const Int_t BLUE[SOL] = { kBlue+2,  kBlue-3,  kBlue-6,  kBlue-9,  kBlue-10, kBlue-10, kBlue-10, kBlue-10,
                                 kBlue-10, kBlue-10, kBlue-10, kBlue-10, kBlue-10, kBlue-10, kBlue-10, kBlue-10 };
static const Int_t GREY[SOL] = { kBlack, kGray+2, kGray+1, kGray,   kGray,   kGray,   kGray,   kGray,
                                 kGray,  kGray,   kGray,   kGray,   kGray,   kGray,   kGray,   kGray };

//-----------------------------------------------------------------------------

enum{ MODEL = 1, BORN = 2, NONRES = 3 };

//-----------------------------------------------------------------------------

void Extract(Char_t* REACT, Int_t L_MAX, Int_t SOLUTIONS=1, Double_t MASS_INITIAL=938.2720)
{
  Double_t MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;

  Char_t Multi[4*(LBINS-1)][256];
  Char_t Fancy[4*(LBINS-1)][256];
  Char_t Buffer1[256];
  Char_t Buffer2[256];
  Double_t Re, DRe, Im, DIm, En;
  Int_t Index;
  TFile* ReFile;
  TFile* ImFile;
  TGraphErrors* ReGraph;
  TGraphErrors* ImGraph;
  FILE* Out;

  //Prepare names for s,p multipoles
  sprintf(Multi[0], "E0p"); sprintf(Fancy[0], "E0+");
  sprintf(Multi[1], "E1p"); sprintf(Fancy[1], "E1+");
  sprintf(Multi[2], "M1p"); sprintf(Fancy[2], "M1+");
  sprintf(Multi[3], "M1m"); sprintf(Fancy[3], "M1-");
  //Prepare names for higher multipoles
  for(Int_t l=2; l<L_MAX+1; l++)
  {
    Index = 4*(l-1);
    sprintf(Multi[Index+0], "E%1dp", l); sprintf(Fancy[Index+0], "E%1d+", l);
    sprintf(Multi[Index+1], "E%1dm", l); sprintf(Fancy[Index+1], "E%1d-", l);
    sprintf(Multi[Index+2], "M%1dp", l); sprintf(Fancy[Index+2], "M%1d+", l);
    sprintf(Multi[Index+3], "M%1dm", l); sprintf(Fancy[Index+3], "M%1d-", l);
  }

  //Process all stored solutions
  for(Int_t s=0; s<SOLUTIONS; s++)
  {
    //Four multipoles (E+, E-, M+, M-) per angular momentum
    for(Int_t t=0; t<4*L_MAX; t++)
    {
      //Open files with real and imaginary part graphs
      sprintf(Buffer1, "plots.%d/Re%s.root", s, Multi[t]);
      ReFile = new TFile(Buffer1);
      sprintf(Buffer1, "plots.%d/Im%s.root", s, Multi[t]);
      ImFile = new TFile(Buffer1);
      //Get graphs for real and imaginary part from files
      sprintf(Buffer1, "Re%s", Multi[t]);
      sprintf(Buffer2, "Re%s_fit", Multi[t]);
      ReGraph = (TGraphErrors*)((TCanvas*)ReFile->Get(Buffer1))->GetPrimitive(Buffer2);
      sprintf(Buffer1, "Im%s", Multi[t]);
      sprintf(Buffer2, "Im%s_fit", Multi[t]);
      ImGraph = (TGraphErrors*)((TCanvas*)ImFile->Get(Buffer1))->GetPrimitive(Buffer2);

      //Prepare output file
      sprintf(Buffer1, "plots.%d/%s.txt", s, Multi[t]);
      Out = fopen(Buffer1, "w");
      //Print header
      fprintf(Out, "   W                 %s(%s)\n", Fancy[t], REACT);
      fprintf(Out, " (MeV)      Re       DRe      Im       DIm\n");
      //Print multipole value for each energy
      for(Int_t n=0; n<ReGraph->GetN(); n++)
      {
        ReGraph->GetPoint(n, En, Re);
        ImGraph->GetPoint(n, En, Im);
        DRe = ReGraph->GetErrorY(n);
        DIm = ImGraph->GetErrorY(n);
        fprintf(Out, "%8.3f  %7.3f  %7.3f  %7.3f  %7.3f\n", TMath::Sqrt(2.0*En*MASS_INITIAL + MASS2_INITIAL), Re, DRe, Im, DIm);
      }
      //Close files
      ReFile->Close();
      ImFile->Close();
      fclose(Out);
    }
  }
}

//-----------------------------------------------------------------------------

void Multipole(Char_t* Mlp, Int_t SolLo=0, Int_t SolHi=0, Bool_t W=true, Bool_t SAVE=false, Double_t Lo=0.0, Double_t Hi=0.0, Double_t MASS_INITIAL=938.2720, Int_t D_WAVES=MODEL)
{
  Double_t MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;

  if(SAVE) TCanvas* Canvas = new TCanvas();
  //if(SAVE) SetBit(TH1::kNoTitle);

  FILE* InPlots;
  FILE* InModel;
  Char_t Buffer[256];
  Char_t M, p;
  Int_t l;
  Double_t PloRe[SOL][N_MAX];
  Double_t PloIm[SOL][N_MAX];
  Double_t PloDRe[SOL][N_MAX];
  Double_t PloDIm[SOL][N_MAX];
  Double_t PloEn[SOL][N_MAX];
  Double_t ModRe[N_MAX];
  Double_t ModIm[N_MAX];
  Double_t ModEn[N_MAX];
  Int_t PloPts[SOL];
  Int_t ModPts = 0;
  Double_t En, Re, DRe, Im, DIm;
  Double_t Min = 0.0;
  Double_t Max = 0.0;
  TGraphErrors* PlotRe[SOL];
  TGraphErrors* PlotIm[SOL];
  TGraph* ModelRe;
  TGraph* ModelIm;

  for(Int_t s=SolLo; s<SolHi+1; s++)
  {
    PloPts[s] = 0;

    //Open text file with fit results for given multipole
    sprintf(Buffer, "plots.%d/%s.txt", s, Mlp);
    InPlots = fopen(Buffer, "r");

    //Skip two lines with table header
    fgets(Buffer, sizeof(Buffer), InPlots);
    fgets(Buffer, sizeof(Buffer), InPlots);
    //Read multipole fit values from file
    while(!feof(InPlots))
    {
      if(fscanf(InPlots, "%lf %lf %lf %lf %lf", &En, &Re, &DRe, &Im, &DIm)==5)
      {
        //Look for maximum & minimum of values to adjust plotting range
        if(Re+DRe > Max) Max = Re+DRe;
        if(Im+DIm > Max) Max = Im+DIm;
        if(Re-DRe < Min) Min = Re-DRe;
        if(Im-DIm < Min) Min = Im-DIm;
        //Plot against either center-of-mass energy W or beam energy omega
        if(!W) En = (En*En - MASS2_INITIAL)/(2.0*MASS_INITIAL);
        //Add real and imaginary parts of multipole (w/ errors) to graph
        PloEn[s][PloPts[s]] = En;
        PloRe[s][PloPts[s]] = Re;
        PloIm[s][PloPts[s]] = Im;
        PloDRe[s][PloPts[s]] = DRe;
        PloDIm[s][PloPts[s]] = DIm;
        PloPts[s]++;
      }
    }
    //Close file with fit values
    fclose(InPlots);

    //Create graphs for real and imaginary parts of fitted multipole
    PlotRe[s] = new TGraphErrors(PloPts[s], PloEn[s], PloRe[s], NULL, PloDRe[s]);
    PlotIm[s] = new TGraphErrors(PloPts[s], PloEn[s], PloIm[s], NULL, PloDIm[s]);

    //Color, line size, marker style adjustments
    PlotRe[s]->SetLineColor(RED[s]);
    PlotIm[s]->SetLineColor(BLUE[s]);
    PlotRe[s]->SetMarkerColor(RED[s]);
    PlotIm[s]->SetMarkerColor(BLUE[s]);
    PlotRe[s]->SetMarkerStyle(21);
    PlotIm[s]->SetMarkerStyle(21);
    PlotRe[s]->SetMarkerSize(0.7);
    PlotIm[s]->SetMarkerSize(0.7);

    //Set plot titles and object names
    PlotRe[s]->SetTitle(Mlp);
    PlotIm[s]->SetTitle(Mlp);
    sprintf(Buffer, "Re%s_fit_%d", Mlp, s);
    PlotRe[s]->SetName(Buffer);
    sprintf(Buffer, "Im%s_fit_%d", Mlp, s);
    PlotIm[s]->SetName(Buffer);

    //Plot graphs
    if(s==SolLo)
    {
      //Adjust x-axis labelling according to W or omega mode
      if(!W)
        PlotRe[s]->GetXaxis()->SetTitle("#omega / MeV");
      else
        PlotRe[s]->GetXaxis()->SetTitle("W / MeV");
      PlotRe[s]->Draw("APZ"); //Plot with x,y axis, points and small error bars
    }
    else
      PlotRe[s]->Draw("PZ"); //Plot with points and small error bars, into existing frame
    PlotIm[s]->Draw("PZ");  //Plot with points and small error bars, into existing frame
  }

  //Open text file with model values for given multipole
  sprintf(Buffer, "model/%s.txt", Mlp);
  sscanf(Mlp, "%c%d%c", &M, &l, &p);
  if((l==2) && (D_WAVES==BORN))   sprintf(Buffer, "model/%s_Born.txt", Mlp);
  if((l==2) && (D_WAVES==NONRES)) sprintf(Buffer, "model/%s_BornRhoOmega.txt", Mlp);
  InModel = fopen(Buffer, "r");

  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), InModel);
  fgets(Buffer, sizeof(Buffer), InModel);
  //Read multipole model values from file
  while(!feof(InModel))
  {
    if(fscanf(InModel, "%lf %lf %lf", &En, &Re, &Im)==3)
    {
      //Look for maximum & minimum of values to adjust plotting range
      if(Re > Max) Max = Re;
      if(Im > Max) Max = Im;
      if(Re < Min) Min = Re;
      if(Im < Min) Min = Im;
      //Plot against either center-of-mass energy W or beam energy omega
      if(!W) En = (En*En - MASS2_INITIAL)/(2.0*MASS_INITIAL);
      //Add real and imaginary parts of multipole to graph
      ModEn[ModPts] = En;
      ModRe[ModPts] = Re;
      ModIm[ModPts] = Im;
      ModPts++;
    }
  }
  //Close file with model values
  fclose(InModel);

  //Create graphs for real and imaginary parts of model multipole
  ModelRe = new TGraph(ModPts, ModEn, ModRe);
  ModelIm = new TGraph(ModPts, ModEn, ModIm);

  //Color, line size adjustments
  ModelRe->SetLineColor(kRed);
  ModelIm->SetLineColor(kBlue);
  ModelRe->SetLineWidth(2);
  ModelIm->SetLineWidth(2);

  //Set plot titles and object names
  ModelRe->SetTitle(Mlp);
  ModelIm->SetTitle(Mlp);
  sprintf(Buffer, "Re%s_model", Mlp);
  ModelRe->SetName(Buffer);
  sprintf(Buffer, "Im%s_model", Mlp);
  ModelIm->SetName(Buffer);

  //Plot graphs
  ModelRe->Draw("L"); //Plot as line, into same frame
  ModelIm->Draw("L"); //Plot as line, into same frame

  //Adjust drawing ranges for y-axis
  if((Lo==0.0) && (Hi==0.0))
  {
    PlotRe[SolLo]->SetMinimum(Min*1.05);
    PlotRe[SolLo]->SetMaximum(Max*1.05);
  }
  else
  {
    PlotRe[SolLo]->SetMinimum(Lo);
    PlotRe[SolLo]->SetMaximum(Hi);
  }

  //y-axis labeling
  sprintf(Buffer, "%c_{%1d%c}", M, l, p);
  for(Int_t n=0; n<strlen(Buffer); n++)
  {
    if(Buffer[n]=='p') Buffer[n] = '+';
    if(Buffer[n]=='m') Buffer[n] = '-';
  }
  PlotRe[SolLo]->GetYaxis()->SetTitle(Buffer);

  if(SAVE) PlotRe[SolLo]->SetTitle("");
  if(SAVE) sprintf(Buffer, "plots.%d/%s.pdf", SolLo, Mlp);
  if(SAVE) Canvas->SaveAs(Buffer);
}

//-----------------------------------------------------------------------------

void Magnitude(Char_t* Mlp, Int_t SolLo=0, Int_t SolHi=0, Bool_t W=true, Bool_t SAVE=false, Double_t Lo=0.0, Double_t Hi=0.0, Double_t MASS_INITIAL=938.2720, Int_t D_WAVES=MODEL)
{
  Double_t MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;

  if(SAVE) TCanvas* Canvas = new TCanvas();
  //if(SAVE) SetBit(TH1::kNoTitle);

  FILE* InPlots;
  FILE* InModel;
  Char_t Buffer[256];
  Char_t M, p;
  Int_t l;
  Double_t PloMg[SOL][N_MAX];
  Double_t PloDMg[SOL][N_MAX];
  Double_t PloEn[SOL][N_MAX];
  Double_t ModMg[N_MAX];
  Double_t ModEn[N_MAX];
  Int_t PloPts[SOL];
  Int_t ModPts = 0;
  Double_t En, Re, DRe, Im, DIm;
  Double_t Min = 0.0;
  Double_t Max = 0.0;
  TGraphErrors* PlotsMg[SOL];
  TGraph* ModelMg;

  for(Int_t s=SolLo; s<SolHi+1; s++)
  {
    PloPts[s] = 0;

    //Open text file with fit results for given multipole
    sprintf(Buffer, "plots.%d/%s.txt", s, Mlp);
    InPlots = fopen(Buffer, "r");

    //Skip two lines with table header
    fgets(Buffer, sizeof(Buffer), InPlots);
    fgets(Buffer, sizeof(Buffer), InPlots);
    //Read multipole fit values from file
    while(!feof(InPlots))
    {
      if(fscanf(InPlots, "%lf %lf %lf %lf %lf", &En, &Re, &DRe, &Im, &DIm)==5)
      {
        //Plot against either center-of-mass energy W or beam energy omega
        if(!W) En = (En*En - MASS2_INITIAL)/(2.0*MASS_INITIAL);
        //Add real and imaginary parts of multipole (w/ errors) to graph
        PloEn[s][PloPts[s]] = En;
        PloMg[s][PloPts[s]] = TMath::Sqrt(Re*Re + Im*Im);
        PloDMg[s][PloPts[s]] = TMath::Sqrt((Re*DRe/PloMg[s][PloPts[s]])*(Re*DRe/PloMg[s][PloPts[s]]) + (Im*DIm/PloMg[s][PloPts[s]])*(Im*DIm/PloMg[s][PloPts[s]]));
        //Look for maximum & minimum of values to adjust plotting range
        if(PloMg[s][PloPts[s]]+PloDMg[s][PloPts[s]] > Max) Max = PloMg[s][PloPts[s]]+PloDMg[s][PloPts[s]];
        if(PloMg[s][PloPts[s]]-PloDMg[s][PloPts[s]] < Min) Min = PloMg[s][PloPts[s]]-PloDMg[s][PloPts[s]];
        PloPts[s]++;
      }
    }
    //Close file with fit values
    fclose(InPlots);

    //Create graphs for real and imaginary parts of fitted multipole
    PlotsMg[s] = new TGraphErrors(PloPts[s], PloEn[s], PloMg[s], NULL, PloDMg[s]);

    //Color, line size, marker style adjustments
    PlotsMg[s]->SetLineColor(GREY[s]);
    PlotsMg[s]->SetMarkerColor(GREY[s]);
    PlotsMg[s]->SetMarkerStyle(21);
    PlotsMg[s]->SetMarkerSize(0.7);

    //Set plot titles and object names
    PlotsMg[s]->SetTitle(Mlp);
    sprintf(Buffer, "Mag%s_fit_%d", Mlp, s);
    PlotsMg[s]->SetName(Buffer);

    //Plot graphs
    if(s==SolLo)
    {
      //Adjust x-axis labelling according to W or omega mode
      if(!W)
        PlotsMg[s]->GetXaxis()->SetTitle("#omega / MeV");
      else
        PlotsMg[s]->GetXaxis()->SetTitle("W / MeV");
      PlotsMg[s]->Draw("APZ"); //Plot with x,y axis, points and small error bars
    }
    else
      PlotsMg[s]->Draw("PZ"); //Plot with points and small error bars, into existing frame
  }

  //Open text file with model values for given multipole
  sprintf(Buffer, "model/%s.txt", Mlp);
  sscanf(Mlp, "%c%d%c", &M, &l, &p);
  if((l==2) && (D_WAVES==BORN))   sprintf(Buffer, "model/%s_Born.txt", Mlp);
  if((l==2) && (D_WAVES==NONRES)) sprintf(Buffer, "model/%s_BornRhoOmega.txt", Mlp);
  InModel = fopen(Buffer, "r");

  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), InModel);
  fgets(Buffer, sizeof(Buffer), InModel);
  //Read multipole model values from file
  while(!feof(InModel))
  {
    if(fscanf(InModel, "%lf %lf %lf", &En, &Re, &Im)==3)
    {
      //Plot against either center-of-mass energy W or beam energy omega
      if(!W) En = (En*En - MASS2_INITIAL)/(2.0*MASS_INITIAL);
      //Add real and imaginary parts of multipole to graph
      ModEn[ModPts] = En;
      ModMg[ModPts] = TMath::Sqrt(Re*Re + Im*Im);
      //Look for maximum & minimum of values to adjust plotting range
      if(ModMg[ModPts] > Max) Max = ModMg[ModPts];
      if(ModMg[ModPts] < Min) Min = ModMg[ModPts];
      ModPts++;
    }
  }
  //Close file with model values
  fclose(InModel);

  //Create graphs for real and imaginary parts of model multipole
  ModelMg = new TGraph(ModPts, ModEn, ModMg);

  //Color, line size adjustments
  ModelMg->SetLineColor(kGray+2);
  ModelMg->SetLineWidth(2);

  //Set plot titles and object names
  ModelMg->SetTitle(Mlp);
  sprintf(Buffer, "Mag%s_model", Mlp);
  ModelMg->SetName(Buffer);

  //Plot graphs
  ModelMg->Draw("L"); //Plot as line, into same frame

  //Adjust drawing ranges for y-axis
  if((Lo==0.0) && (Hi==0.0))
  {
    PlotsMg[0]->SetMinimum(Min*1.05);
    PlotsMg[0]->SetMaximum(Max*1.05);
  }
  else
  {
    PlotsMg[0]->SetMinimum(Lo);
    PlotsMg[0]->SetMaximum(Hi);
  }

  //y-axis labeling
  sprintf(Buffer, "|%c_{%1d%c}|", M, l, p);
  for(Int_t n=0; n<strlen(Buffer); n++)
  {
    if(Buffer[n]=='p') Buffer[n] = '+';
    if(Buffer[n]=='m') Buffer[n] = '-';
  }
  PlotsMg[0]->GetYaxis()->SetTitle(Buffer);

  if(SAVE) PlotsMg[0]->SetTitle("");
  if(SAVE) sprintf(Buffer, "plots.0/Mag%s.pdf", Mlp);
  if(SAVE) Canvas->SaveAs(Buffer);
}
//-----------------------------------------------------------------------------

void Phase(Char_t* Mlp, Int_t SolLo=0, Int_t SolHi=0, Bool_t W=true, Bool_t SAVE=false, Double_t Lo=0.0, Double_t Hi=0.0, Double_t MASS_INITIAL=938.2720, Int_t D_WAVES=MODEL)
{
  Double_t MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;

  if(SAVE) TCanvas* Canvas = new TCanvas();
  //if(SAVE) SetBit(TH1::kNoTitle);

  FILE* InPlots;
  FILE* InModel;
  Char_t Buffer[256];
  Char_t M, p;
  Int_t l;
  Double_t PloPh[SOL][N_MAX];
  Double_t PloDPh[SOL][N_MAX];
  Double_t PloEn[SOL][N_MAX];
  Double_t ModPh[N_MAX];
  Double_t ModEn[N_MAX];
  Int_t PloPts[SOL];
  Int_t ModPts = 0;
  Double_t En, Re, DRe, Im, DIm;
  Double_t dz;
  Double_t Min = 0.0;
  Double_t Max = 0.0;
  TGraphErrors* PlotsPh[SOL];
  TGraph* ModelPh;

  for(Int_t s=SolLo; s<SolHi+1; s++)
  {
    PloPts[s] = 0;

    //Open text file with fit results for given multipole
    sprintf(Buffer, "plots.%d/%s.txt", s, Mlp);
    InPlots = fopen(Buffer, "r");

    //Skip two lines with table header
    fgets(Buffer, sizeof(Buffer), InPlots);
    fgets(Buffer, sizeof(Buffer), InPlots);
    //Read multipole fit values from file
    while(!feof(InPlots))
    {
      if(fscanf(InPlots, "%lf %lf %lf %lf %lf", &En, &Re, &DRe, &Im, &DIm)==5)
      {
        //Plot against either center-of-mass energy W or beam energy omega
        if(!W) En = (En*En - MASS2_INITIAL)/(2.0*MASS_INITIAL);
        //Add real and imaginary parts of multipole (w/ errors) to graph
        PloEn[s][PloPts[s]] = En;
        PloPh[s][PloPts[s]] = TMath::ATan(Im/Re);
        dz = 1.0/(1.0 + (Im*Im)/(Re*Re));
        PloDPh[s][PloPts[s]] = TMath::Sqrt((dz*DIm/Re)*(dz*DIm/Re) + (dz*DRe*Im/(Re*Re))*(dz*DRe*Im/(Re*Re)));
        //Look for maximum & minimum of values to adjust plotting range
        if(PloPh[s][PloPts[s]]+PloDPh[s][PloPts[s]] > Max) Max = PloPh[s][PloPts[s]]+PloDPh[s][PloPts[s]];
        if(PloPh[s][PloPts[s]]-PloDPh[s][PloPts[s]] < Min) Min = PloPh[s][PloPts[s]]-PloDPh[s][PloPts[s]];
        PloPts[s]++;
      }
    }
    //Close file with fit values
    fclose(InPlots);

    //Create graphs for real and imaginary parts of fitted multipole
    PlotsPh[s] = new TGraphErrors(PloPts[s], PloEn[s], PloPh[s], NULL, PloDPh[s]);

    //Color, line size, marker style adjustments
    PlotsPh[s]->SetLineColor(GREY[s]);
    PlotsPh[s]->SetMarkerColor(GREY[s]);
    PlotsPh[s]->SetMarkerStyle(21);
    PlotsPh[s]->SetMarkerSize(0.7);

    //Set plot titles and object names
    PlotsPh[s]->SetTitle(Mlp);
    sprintf(Buffer, "Phs%s_fit_%d", Mlp, s);
    PlotsPh[s]->SetName(Buffer);

    //Plot graphs
    if(s==0)
    {
      //Adjust x-axis labelling according to W or omega mode
      if(!W)
        PlotsPh[s]->GetXaxis()->SetTitle("#omega / MeV");
      else
        PlotsPh[s]->GetXaxis()->SetTitle("W / MeV");
      PlotsPh[s]->Draw("APZ"); //Plot with x,y axis, points and small error bars
    }
    else
      PlotsPh[s]->Draw("PZ"); //Plot with points and small error bars, into existing frame
  }

  //Open text file with model values for given multipole
  sprintf(Buffer, "model/%s.txt", Mlp);
  sscanf(Mlp, "%c%d%c", &M, &l, &p);
  if((l==2) && (D_WAVES==BORN))   sprintf(Buffer, "model/%s_Born.txt", Mlp);
  if((l==2) && (D_WAVES==NONRES)) sprintf(Buffer, "model/%s_BornRhoOmega.txt", Mlp);
  InModel = fopen(Buffer, "r");

  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), InModel);
  fgets(Buffer, sizeof(Buffer), InModel);
  //Read multipole model values from file
  while(!feof(InModel))
  {
    if(fscanf(InModel, "%lf %lf %lf", &En, &Re, &Im)==3)
    {
      //Plot against either center-of-mass energy W or beam energy omega
      if(!W) En = (En*En - MASS2_INITIAL)/(2.0*MASS_INITIAL);
      //Add real and imaginary parts of multipole to graph
      ModEn[ModPts] = En;
      ModPh[ModPts] = TMath::ATan(Im/Re);
      //Look for maximum & minimum of values to adjust plotting range
      if(ModPh[ModPts] > Max) Max = ModPh[ModPts];
      if(ModPh[ModPts] < Min) Min = ModPh[ModPts];
      ModPts++;
    }
  }
  //Close file with model values
  fclose(InModel);

  //Create graphs for real and imaginary parts of model multipole
  ModelPh = new TGraph(ModPts, ModEn, ModPh);

  //Color, line size adjustments
  ModelPh->SetLineColor(kGray+2);
  ModelPh->SetLineWidth(2);

  //Set plot titles and object names
  ModelPh->SetTitle(Mlp);
  sprintf(Buffer, "Phs%s_model", Mlp);
  ModelPh->SetName(Buffer);

  //Plot graphs
  ModelPh->Draw("L"); //Plot as line, into same frame

  //Adjust drawing ranges for y-axis
  if((Lo==0.0) && (Hi==0.0))
  {
    PlotsPh[0]->SetMinimum(Min*1.05);
    PlotsPh[0]->SetMaximum(Max*1.05);
  }
  else
  {
    PlotsPh[0]->SetMinimum(Lo);
    PlotsPh[0]->SetMaximum(Hi);
  }

  //y-axis labeling
  sprintf(Buffer, "#phi_{%c_{%1d%c}}", M, l, p);
  for(Int_t n=4; n<strlen(Buffer); n++)
  {
    if(Buffer[n]=='p') Buffer[n] = '+';
    if(Buffer[n]=='m') Buffer[n] = '-';
  }
  PlotsPh[0]->GetYaxis()->SetTitle(Buffer);

  if(SAVE) PlotsPh[0]->SetTitle("");
  if(SAVE) sprintf(Buffer, "plots.0/Phs%s.pdf", Mlp);
  if(SAVE) Canvas->SaveAs(Buffer);
}

//-----------------------------------------------------------------------------

void Model(Char_t* Mlp, Bool_t W=true, Bool_t SAVE=false, Double_t Lo=0.0, Double_t Hi=0.0, Double_t MASS_INITIAL=938.2720, Int_t D_WAVES=MODEL)
{
  Double_t MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;

  if(SAVE) TCanvas* Canvas = new TCanvas();
  //if(SAVE) SetBit(TH1::kNoTitle);

  FILE* InModel;
  Char_t Buffer[256];
  Char_t M, p;
  Int_t l;
  Double_t ModRe[N_MAX];
  Double_t ModIm[N_MAX];
  Double_t ModEn[N_MAX];
  Int_t ModPts = 0;
  Double_t En, Re, Im;
  Double_t Min = 0.0;
  Double_t Max = 0.0;
  TGraph* ModelRe;
  TGraph* ModelIm;

  //Open text file with model values for given multipole
  sprintf(Buffer, "model/%s.txt", Mlp);
  sscanf(Mlp, "%c%d%c", &M, &l, &p);
  if((l==2) && (D_WAVES==BORN))   sprintf(Buffer, "model/%s_Born.txt", Mlp);
  if((l==2) && (D_WAVES==NONRES)) sprintf(Buffer, "model/%s_BornRhoOmega.txt", Mlp);
  InModel = fopen(Buffer, "r");

  //Skip two lines with table header
  fgets(Buffer, sizeof(Buffer), InModel);
  fgets(Buffer, sizeof(Buffer), InModel);
  //Read multipole model values from file
  while(!feof(InModel))
  {
    if(fscanf(InModel, "%lf %lf %lf", &En, &Re, &Im)==3)
    {
      //Look for maximum & minimum of values to adjust plotting range
      if(Re > Max) Max = Re;
      if(Im > Max) Max = Im;
      if(Re < Min) Min = Re;
      if(Im < Min) Min = Im;
      //Plot against either center-of-mass energy W or beam energy omega
      if(!W) En = (En*En - MASS2_INITIAL)/(2.0*MASS_INITIAL);
      //Add real and imaginary parts of multipole to graph
      ModEn[ModPts] = En;
      ModRe[ModPts] = Re;
      ModIm[ModPts] = Im;
      ModPts++;
    }
  }
  //Close file with model values
  fclose(InModel);

  //Create graphs for real and imaginary parts of model multipole
  ModelRe = new TGraph(ModPts, ModEn, ModRe);
  ModelIm = new TGraph(ModPts, ModEn, ModIm);

  //Color, line size adjustments
  ModelRe->SetLineColor(kRed);
  ModelIm->SetLineColor(kBlue);
  ModelRe->SetLineWidth(2);
  ModelIm->SetLineWidth(2);

  //Set plot titles and object names
  ModelRe->SetTitle(Mlp);
  ModelIm->SetTitle(Mlp);
  sprintf(Buffer, "Re%s_model", Mlp);
  ModelRe->SetName(Buffer);
  sprintf(Buffer, "Im%s_model", Mlp);
  ModelIm->SetName(Buffer);

  //Plot graphs
  ModelRe->Draw("L"); //Plot as line, into same frame
  ModelIm->Draw("L"); //Plot as line, into same frame

  //Adjust drawing ranges for y-axis
  if((Lo==0.0) && (Hi==0.0))
  {
    ModelRe->SetMinimum(Min*1.05);
    ModelRe->SetMaximum(Max*1.05);
  }
  else
  {
    ModelRe->SetMinimum(Lo);
    ModelRe->SetMaximum(Hi);
  }

  //y-axis labeling
  sprintf(Buffer, "%c_{%1d%c}", M, l, p);
  for(Int_t n=0; n<strlen(Buffer); n++)
  {
    if(Buffer[n]=='p') Buffer[n] = '+';
    if(Buffer[n]=='m') Buffer[n] = '-';
  }
  ModelRe->GetYaxis()->SetTitle(Buffer);

  //Plot graphs
  ModelRe->Draw("AL"); //Plot as line, with x,y axis
  ModelIm->Draw("L");  //Plot with line into existing frame
  //Adjust x-axis labelling according to W or omega mode
  if(!W)
    ModelRe->GetXaxis()->SetTitle("#omega / MeV");
  else
    ModelRe->GetXaxis()->SetTitle("W / MeV");

  if(SAVE) ModelRe->SetTitle("");
  if(SAVE) sprintf(Buffer, "%s.pdf", Mlp);
  if(SAVE) Canvas->SaveAs(Buffer);
}

//-----------------------------------------------------------------------------

void Chi2(Int_t SOLUTION=0, Double_t MASS_INITIAL=938.2720)
{
  Double_t MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;

  Char_t Buffer[256];
  Double_t InE, InChi2;
  TFile* Chi2File;
  TGraph* Chi2GraphE;
  TGraph* Chi2GraphW;
  TCanvas* Canvas;
  Double_t Chi2[N_MAX];
  Double_t W[N_MAX];
  Int_t Chi2Pts;

  sprintf(Buffer, "plots.%d/Chi2.root", SOLUTION);
  Chi2File = new TFile(Buffer);
  Chi2GraphE = (TGraph*)((TCanvas*)Chi2File->Get("Chi2"))->GetPrimitive("Chi2");

  Chi2Pts = Chi2GraphE->GetN();
  for(Int_t n=0; n<Chi2Pts; n++)
  {
    Chi2GraphE->GetPoint(n, InE, InChi2);
    W[n] = TMath::Sqrt(2.0*InE*MASS_INITIAL + MASS2_INITIAL);
    Chi2[n] = InChi2;
  }
  Chi2File->Close();

  Canvas = new TCanvas();
  Chi2GraphW = new TGraph(Chi2Pts, W, Chi2);
  Chi2GraphW->SetMarkerStyle(21);
  Chi2GraphW->SetMarkerSize(0.7);
  Chi2GraphW->SetName("Chi2");
  Chi2GraphW->SetTitle("Chi2");
  Chi2GraphW->Draw("APZ");
}

//-----------------------------------------------------------------------------

void Penalty(Int_t SOLUTION=0, Double_t MASS_INITIAL=938.2720)
{
  Double_t MASS2_INITIAL = MASS_INITIAL*MASS_INITIAL;

  Char_t Buffer[256];
  Double_t InE, InPen;
  TFile* PenFile;
  TGraph* PenGraphE;
  TGraph* PenGraphW;
  TCanvas* Canvas;
  Double_t Pen[N_MAX];
  Double_t W[N_MAX];
  Int_t PenPts;

  sprintf(Buffer, "plots.%d/Penalty.root", SOLUTION);
  PenFile = new TFile(Buffer);
  PenGraphE = (TGraph*)((TCanvas*)PenFile->Get("Penalty"))->GetPrimitive("Penalty");

  PenPts = PenGraphE->GetN();
  for(Int_t n=0; n<PenPts; n++)
  {
    PenGraphE->GetPoint(n, InE, InPen);
    W[n] = TMath::Sqrt(2.0*InE*MASS_INITIAL + MASS2_INITIAL);
    Pen[n] = InPen;
  }
  PenFile->Close();

  Canvas = new TCanvas();
  PenGraphW = new TGraph(PenPts, W, Pen);
  PenGraphW->SetMarkerStyle(21);
  PenGraphW->SetMarkerSize(0.7);
  PenGraphW->SetName("Penalty");
  PenGraphW->SetTitle("Penalty");
  PenGraphW->Draw("APZ");
}

//-----------------------------------------------------------------------------
