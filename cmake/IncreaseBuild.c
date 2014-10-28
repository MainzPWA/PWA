#include <stdio.h>

int main( int argc, const char* argv[] )
{
  int build = 0;
  FILE* file = NULL;

  file = fopen(argv[1], "r");
  if(file)
  {
    fscanf(file, "#define BUILD %d\n", &build);
    fclose(file);
  }

  file = fopen(argv[1], "w");
  if(file)
  {
    fprintf(file, "#define BUILD %d\n", build+1);
    fclose(file);
  }
}
