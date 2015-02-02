void ChMAIDConv(Char_t* Mlp)
{
  Char_t Buffer[256];
  TMultiDim fint_re(1);
  TMultiDim fint_im(1);
  Double_t En, Re, Im;
  Int_t Idx = 0;
  FILE* In;
  FILE* Out;

  fint_re.SetNBins(128);
  fint_im.SetNBins(128);

  sprintf(Buffer, "%s.txt", Mlp);
  In = fopen(Buffer, "r");
  while(!feof(In))
  {
    fscanf(In, "%lf %lf %lf\n", &En, &Re, &Im);
    fint_re.SetValue(Idx, Re);
    fint_im.SetValue(Idx, Im);
    fint_re.SetCoord(0, Idx, En);
    fint_im.SetCoord(0, Idx, En);
    Idx++;
  }
  fclose(In);

  sprintf(Buffer, "%s_new.txt", Mlp);
  Out = fopen(Buffer, "w");
  for(Int_t t=0; t<strlen(Mlp); t++)
  {
    if(Mlp[t]=='p') Mlp[t] = '+';
    if(Mlp[t]=='m') Mlp[t] = '-';
  }
  fprintf(Out, "   W          %s(pi0_p)\n", Mlp);
  fprintf(Out, " (MeV)       Re         Im\n");
  for(Double_t e=145.0; e<200.1; e+=0.5)
  {
    En = TMath::Sqrt(2.0*938.272*e + 938.272*938.272);
    fprintf(Out, "%8.3f  %8.4f  %8.4f\n",
            En, fint_re.Interpolate(En), fint_im.Interpolate(En));
  }
  fclose(Out);
}
