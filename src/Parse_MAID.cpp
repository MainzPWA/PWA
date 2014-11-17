#include "Parse_MAID.h"

//-----------------------------------------------------------------------------

void Parse_MAID(Char_t* Path)
{
  Char_t Buffer[1024];
  Double_t Energy, Re, Im;
  FILE* MAID_Ep[LBINS];
  FILE* MAID_Em[LBINS];
  FILE* MAID_Mp[LBINS];
  FILE* MAID_Mm[LBINS];

  EPSILON = 1e38; //Initialize for search

  printf("------------------------------------------------------------------------------------\n");
  printf("Loading model multipoles... ");

  //Open s,p waves
  sprintf(Buffer, "%s/E0p.txt", Path); MAID_Ep[0] = fopen(Buffer, "r");
  sprintf(Buffer, "%s/E1p.txt", Path); MAID_Ep[1] = fopen(Buffer, "r");
  sprintf(Buffer, "%s/M1p.txt", Path); MAID_Mp[1] = fopen(Buffer, "r");
  sprintf(Buffer, "%s/M1m.txt", Path); MAID_Mm[1] = fopen(Buffer, "r");
  //Open d waves (depending on parametrisation)
  switch(D_WAVES)
  {
   case BORN: //Born terms only (in the strict sense)
    sprintf(Buffer, "%s/E2p_Born.txt", Path); MAID_Ep[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/E2m_Born.txt", Path); MAID_Em[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/M2p_Born.txt", Path); MAID_Mp[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/M2m_Born.txt", Path); MAID_Mm[2] = fopen(Buffer, "r");
    break;
   case NONRES: //Born, rho, omega terms (=non-resonant contributions)
    sprintf(Buffer, "%s/E2p_BornRhoOmega.txt", Path); MAID_Ep[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/E2m_BornRhoOmega.txt", Path); MAID_Em[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/M2p_BornRhoOmega.txt", Path); MAID_Mp[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/M2m_BornRhoOmega.txt", Path); MAID_Mm[2] = fopen(Buffer, "r");
    break;
   default: //Full model prediction
    sprintf(Buffer, "%s/E2p.txt", Path); MAID_Ep[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/E2m.txt", Path); MAID_Em[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/M2p.txt", Path); MAID_Mp[2] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/M2m.txt", Path); MAID_Mm[2] = fopen(Buffer, "r");
    break;
  }
  //Open f,g,... waves
  for(Int_t l=3; l<L_MAX+1; l++)
  {
    sprintf(Buffer, "%s/E%dp.txt", Path, l); MAID_Ep[l] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/E%dm.txt", Path, l); MAID_Em[l] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/M%dp.txt", Path, l); MAID_Mp[l] = fopen(Buffer, "r");
    sprintf(Buffer, "%s/M%dm.txt", Path, l); MAID_Mm[l] = fopen(Buffer, "r");
  }

  //Load E0+ multipole
  maid_bin = 0;
  //Skip 2 uninteresting lines
  fgets(Buffer, sizeof(Buffer), MAID_Ep[0]);
  fgets(Buffer, sizeof(Buffer), MAID_Ep[0]);
  while(!feof(MAID_Ep[0]))
  {
    if(fscanf(MAID_Ep[0], "%lf %lf %lf\n", &Energy, &Re, &Im)!=3) break;
    //Find absolute non-zero minimum of Re or Im in model prediction
    if((Abs(Re) < EPSILON) && (Abs(Re) > 0.0)) EPSILON = Abs(Re);
    if((Abs(Im) < EPSILON) && (Abs(Im) > 0.0)) EPSILON = Abs(Im);
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
    //Find absolute non-zero minimum of Re or Im in model prediction
    if((Abs(Re) < EPSILON) && (Abs(Re) > 0.0)) EPSILON = Abs(Re);
    if((Abs(Im) < EPSILON) && (Abs(Im) > 0.0)) EPSILON = Abs(Im);
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
    //Find absolute non-zero minimum of Re or Im in model prediction
    if((Abs(Re) < EPSILON) && (Abs(Re) > 0.0)) EPSILON = Abs(Re);
    if((Abs(Im) < EPSILON) && (Abs(Im) > 0.0)) EPSILON = Abs(Im);
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
    //Find absolute non-zero minimum of Re or Im in model prediction
    if((Abs(Re) < EPSILON) && (Abs(Re) > 0.0)) EPSILON = Abs(Re);
    if((Abs(Im) < EPSILON) && (Abs(Im) > 0.0)) EPSILON = Abs(Im);
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
      //Find absolute non-zero minimum of Re or Im in model prediction
      if((Abs(Re) < EPSILON) && (Abs(Re) > 0.0)) EPSILON = Abs(Re);
      if((Abs(Im) < EPSILON) && (Abs(Im) > 0.0)) EPSILON = Abs(Im);
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
      //Find absolute non-zero minimum of Re or Im in model prediction
      if((Abs(Re) < EPSILON) && (Abs(Re) > 0.0)) EPSILON = Abs(Re);
      if((Abs(Im) < EPSILON) && (Abs(Im) > 0.0)) EPSILON = Abs(Im);
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
      //Find absolute non-zero minimum of Re or Im in model prediction
      if((Abs(Re) < EPSILON) && (Abs(Re) > 0.0)) EPSILON = Abs(Re);
      if((Abs(Im) < EPSILON) && (Abs(Im) > 0.0)) EPSILON = Abs(Im);
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
      //Find absolute non-zero minimum of Re or Im in model prediction
      if((Abs(Re) < EPSILON) && (Abs(Re) > 0.0)) EPSILON = Abs(Re);
      if((Abs(Im) < EPSILON) && (Abs(Im) > 0.0)) EPSILON = Abs(Im);
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

  //Make EPSILON a small number, i.e. 50% of the found absolute non-zero minimum of
  //model predictions, corresponding to half a digit of the model's precision.
  EPSILON*=0.5;
  EPSILON2 = EPSILON*EPSILON;

  printf("%2d multipoles at %4d energies loaded\n", 4*L_MAX, maid_bin);
  printf("------------------------------------------------------------------------------------\n");
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

Int_t GetEnergyBin_maid()
{
  //Get energy bin for sigma0 for given global energy
  Double_t Min = 1e38;
  Int_t eM = 0;

  for(Int_t e=0; e<maid_bin; e++)
    if(fabs(maid_en[e] - gEnergy) < Min)
    {
      Min = fabs(maid_en[e] - gEnergy);
      eM = e;
    }
  return eM;
}

//-----------------------------------------------------------------------------

