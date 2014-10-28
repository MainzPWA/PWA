static const Int_t EBINS  = 2048;
static const Int_t THBINS = 256;

//-----------------------------------------------------------------------------

Double_t Obs_val[EBINS][THBINS];
Double_t Obs_err[EBINS][THBINS];
Double_t Obs_unc[EBINS][THBINS];
Double_t Obs_th[EBINS][THBINS];
Double_t Obs_en[EBINS];
Double_t Obs_wt[EBINS];
Double_t Obs_sy[EBINS];
Int_t Obs_pts[EBINS];
Int_t Obs_bin;

//-----------------------------------------------------------------------------

void Parse_Obs(Char_t* InFile)
{
  Char_t Buffer[1024];
  Int_t ThetaBin, EnergyBin;
  Double_t Energy, Weight, System, Theta, Obs, DObs, EObs;

  FILE* File_in = fopen(InFile, "r");

  Obs_bin = 0;
  while(!feof(File_in))
  {
    //Get beam energy
    if(fscanf(File_in, "E = %lf MeV, Wght = %lf, Syst = %lf\n", &Energy, &Weight, &System)!=3) break;

    ThetaBin = 0;
    while(fscanf(File_in, "%lf %lf %lf %lf\n", &Theta, &Obs, &DObs, &EObs)==4)
    {
      Obs_val[Obs_bin][ThetaBin] = Obs;
      Obs_err[Obs_bin][ThetaBin] = DObs;
      Obs_unc[Obs_bin][ThetaBin] = EObs;
      Obs_th[Obs_bin][ThetaBin]  = Theta;
      if(DObs!=0.0) ThetaBin++; //Accept only 'existing' data points
    }
    //Skip 1 uninteresting line
    fgets(Buffer, sizeof(Buffer), File_in);

    Obs_pts[Obs_bin] = ThetaBin;
    Obs_en[Obs_bin] = Energy;
    Obs_wt[Obs_bin] = Weight;
    Obs_sy[Obs_bin] = System;

    Obs_bin++;
  }

  fclose(File_in);
  //Count data points and (used) energy bins
  Int_t n = 0; for(Int_t t=0; t<Obs_bin; t++) n+=Obs_pts[t];
  Int_t m = 0; for(Int_t t=0; t<Obs_bin; t++) if(Obs_pts[t]) m++;
  printf("%5d data points at %3d energies loaded\n", n, m);

  /*//Debug output
  printf("EBins: %d\n", Obs_bin);
  for(Int_t e=0; e<Obs_bin; e++)
  {
    printf("%d (%f MeV): ThBins: %d\n", e, Obs_en[e], Obs_pts[e]);
    for(Int_t th=0; th<Obs_pts[e]; th++)
      printf("%f %f %f\n", Obs_th[e][th], Obs_val[e][th], Obs_err[e][th]);
  }*/

  sprintf(Buffer, "out_%s", InFile);
  FILE* File_out = fopen(Buffer, "w");
  for(Int_t e=0; e<Obs_bin; e++)
  {
    fprintf(File_Obs_out, "E = %8.3f MeV, Wght = %4.2f, Syst = %8.6f\n", Obs_en[e], Obs_wt[e], Obs_sy[e]);
    for(Int_t th=0; th<Obs_pts[e]; th++)
      fprintf(File_Obs_out, "%7.3f   %10.6f   %10.6f   %10.6f\n", Obs_th[e][th], Obs_val[e][th], Obs_err[e][th], Obs_unc[e][th]);
    fprintf(File_Obs_out, "----------------------------------------------\n");
  }
  fclose(File_Obs_out);
}

//-----------------------------------------------------------------------------
