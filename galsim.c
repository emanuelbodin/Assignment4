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

typedef struct _particleBox
{
  double xCenter;
  double yCenter;
  double width;
  double length;
  particle star;
} particleBox;

particle * read_particle(int N, FILE *fp1) {
  unsigned char buffer[8];
  double arr[6];
  particle *array = malloc(N * sizeof(particle));
  for(int i = 0; i<N; i++){
    for(int j = 0; j<6; j++){
        fread(buffer, sizeof(buffer), 1, fp1);
        arr[j] = *((double*)buffer);
    }
    array[i].posX = arr[0];
    array[i].posY = arr[1];
    array[i].mass = arr[2];
    array[i].velX = arr[3];
    array[i].velY = arr[4];
    array[i].b = arr[5];
  }

  fclose(fp1);
  return array;
}

// ./galsim 2 input/cicles_N_2.gal 500 1e-5 0.1 0
int main(int argc, char* argv[]){
  if (argc != 7){
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
  
  particle *array = read_particle(N, fp1);

  return 0;
}