void Align(Char_t* Mlp, Char_t* React)
{
  FILE* In;
  FILE* Out;
  Char_t Buffer[256];
  Double_t W, Re, Im;

  sprintf(Buffer, "tmp/%s.txt", Mlp);
  In  = fopen(Buffer, "r");
  sprintf(Buffer, "%s.txt", Mlp);
  Out = fopen(Buffer, "w");
  for(Int_t i=0; i<strlen(Mlp); i++)
  {
    if(Mlp[i]=='p') Mlp[i] = '+';
    if(Mlp[i]=='m') Mlp[i] = '-';
  }
  fprintf(Out, "   W         %s(%s)\n", Mlp, React);
  fprintf(Out, " (MeV)      Re         Im\n");
  while(!feof(In))
  {
    if(fscanf(In, "%lf %lf %lf", &W, &Re, &Im)==3)
    {
      fprintf(Out, " %6.1f  %9.5f  %9.5f\n", W, Re, Im);
    }
  }
  fclose(In);
  fclose(Out);
}
