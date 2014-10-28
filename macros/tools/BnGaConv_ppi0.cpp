void BnGaConv(Char_t* Mlp)
{
  Double_t Conversion = 139.5702/197.3269; //m_pi+/hbar*c converts from fm to 1/m_pi+
  Char_t Buffer[256];
  TMultiDim fint_re(1);
  TMultiDim fint_im(1);
  Double_t En, Re, Im;
  Int_t Idx = 0;
  FILE* In;
  FILE* Out;

  fint_re.SetNBins(512);
  fint_im.SetNBins(512);

  sprintf(Buffer, "pion_00_%s.dat", Mlp);
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

  sprintf(Buffer, "%s.txt", Mlp);
  Out = fopen(Buffer, "w");
  for(Int_t t=0; t<strlen(Mlp); t++)
  {
    if(Mlp[t]=='p') Mlp[t] = '+';
    if(Mlp[t]=='m') Mlp[t] = '-';
  }
  fprintf(Out, "   W         %s(pi0_p)\n", Mlp);
  fprintf(Out, " (MeV)      Re         Im\n");
  for(Double_t e=1074.0; e<1972.0; e+=1.0)
  {
    fprintf(Out, "%7.1f  %10.6f %10.6f\n",
            e, fint_re.Interpolate(e)*Conversion, fint_im.Interpolate(e)*Conversion);
  }
  fclose(Out);
}
