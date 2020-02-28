#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct _particle
{
  double posX;
  double posY;
  double mass;
  double velX;
  double velY;
  double b;
}particle;

// ./galsim 2 input_data/circles_N_2.gagraphics
// gcc -O1 -g -pg galsim.c -lm
int main(int argc, char* argv[]){
  if (argc != 6){
      printf("Wrong number of input arguments\n");
      return 1;
  }
  int N = atoi(argv[1]);
  char* filename = argv[2];
  int n_steps = atoi(argv[3]);
  double delta_t = atof(argv[4]);
  double theta_max = atof(argv[5]);
  int graphics = atoi(argv[6]);

  printf("Command line arguments given: %d, %s, %d, %f, %d, %d \n", N, filename, n_steps, delta_t, theta_max, graphics);
    FILE *fp1, *fp2;
    fp1 = fopen(filename, "rb");
    const double e0 = 0.001;
    const double G = 100.0 / N;

    if (fp1 == NULL){
    printf("Error while opening the file.\n");
    return 1;
    }

  return 0;
}