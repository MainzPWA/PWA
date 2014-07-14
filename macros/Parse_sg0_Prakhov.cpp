static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t sg0_val[EBINS][THBINS];
Double_t sg0_err[EBINS][THBINS];
Double_t sg0_th[EBINS][THBINS];
Double_t sg0_en[EBINS];
Int_t sg0_pts[EBINS];
Int_t sg0_bin;

//-----------------------------------------------------------------------------

void Parse_sg0(Char_t* PrakhovFile, Double_t Weight=1.00, Double_t System=0.05)
{
  Char_t Buffer[1024];
  Int_t ThetaBin;
  Double_t Energy, Theta, sigma0, Dsigma0;
  Double_t Dummy1, Dummy2, Dummy3;
  FILE* File_sg0_in;
  FILE* File_sg0_out;

  File_sg0_in = fopen(PrakhovFile, "r");

  sg0_bin = 0;
  while(!feof(File_sg0_in))
  {
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sg0_in);
    fgets(Buffer, sizeof(Buffer), File_sg0_in);
    //Get beam energy
    if(fscanf(File_sg0_in, "Eg=%lf+/-%lf W=%lf+/-%lf\n", &Energy, &Dummy1, &Dummy2, &Dummy3)!=4) break;
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_sg0_in);
    ThetaBin = 0;
    while(fscanf(File_sg0_in, "%lf %lf %lf %lf %lf\n", &Dummy1, &Theta, &sigma0, &Dsigma0, &Dummy2)==5)
    {
      sg0_val[sg0_bin][ThetaBin] = sigma0;
      sg0_err[sg0_bin][ThetaBin] = Dsigma0;
      sg0_th[sg0_bin][ThetaBin]  = Theta;
      if((sigma0!=0.0) && (Dsigma0!=0.0)) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 2 uninteresting lines
    fgets(Buffer, sizeof(Buffer), File_sg0_in);
    fgets(Buffer, sizeof(Buffer), File_sg0_in);

    sg0_pts[sg0_bin] = ThetaBin;
    sg0_en[sg0_bin] = Energy;
    sg0_bin++;
  }
  fclose(File_sg0_in);

  sprintf(Buffer, "out_%s", PrakhovFile);
  File_sg0_out = fopen(Buffer, "w");
  for(Int_t e=0; e<sg0_bin; e++)
  {
    fprintf(File_sg0_out, "E = %8.3f MeV, Wght = %3.1f, Syst = %6.4f\n", sg0_en[e], Weight, System);
    for(Int_t th=0; th<sg0_pts[e]; th++)
      fprintf(File_sg0_out, "%7.3f  %10.6f  %10.6f\n", sg0_th[e][th], sg0_val[e][th], sg0_err[e][th]);
    fprintf(File_sg0_out, "-------------------------------------------\n");
  }
  fclose(File_sg0_out);
}

//-----------------------------------------------------------------------------
