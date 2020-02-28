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

const float circleRadius=0.0025, circleColor=0;
const int windowWidth=800;

void keep_within_box(double* xA, double* yA) {
  if(*xA > 1)
    *xA = 0;
  if(*yA > 1)
    *yA = 0;
}

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
  printf("%f \n", delta_t);
  int graphics = atoi(argv[5]);

  printf("Command line arguments given: %d, %s, %d, %f, %d \n", N, filename, n_steps, delta_t, graphics);
    FILE *fp1, *fp2;
    fp1 = fopen(filename, "rb");
    const double e0 = 0.001;
    const double G = 100.0 / N;

    if (fp1 == NULL){
    printf("Error while opening the file.\n");
    return 1;
    }
    particle *array = read_particle(N, fp1);

    particle *new_array = malloc(N * sizeof(particle));
    for (int k = 0; k < n_steps; k++) {
      for (int i = 0; i < N; i++) {
        double Fix = 0;
        double Fiy = 0;
        double aix, aiy;
        double uix, uiy;
        double pix, piy;
        double massi = array[i].mass;
        for (int j = 0;  j < N; j++) {
          double Fjx = 0, Fjy = 0;
          if (j == i) {} /* skip iteration */
          else {
            double x_dist = array[i].posX - array[j].posX;
            double y_dist = array[i].posY - array[j].posY;
            double r2 = pow(x_dist, 2) + pow(y_dist, 2);
            double denom = pow((sqrt(r2)+e0), 3);
            double massj = array[j].mass;
            Fjx = massj * x_dist / denom; 
            Fjy = massj * y_dist / denom;
            Fix += Fjx;
            Fiy += Fjy;
          }
        }
        Fix = -G * massi * Fix;
        Fiy = -G * massi * Fiy;
        aix = Fix / massi;
        aiy = Fiy / massi;
        uix = array[i].velX + delta_t * aix;
        uiy = array[i].velY + delta_t * aiy;
        pix = array[i].posX + delta_t * uix;
        piy = array[i].posY + delta_t * uiy;
        // update for next time step
        new_array[i].velX = uix;
        new_array[i].velY = uiy;
        new_array[i].posX = pix;
        new_array[i].posY = piy;
        new_array[i].mass = array[i].mass;
        new_array[i].b = array[i].b;
      }
      
      for(int i = 0; i<N; i++){
        array[i] = new_array[i];
      }
    }
  
    fp2 = fopen("result.gal", "wb");
    for(int i = 0; i<N; i++){
      fwrite(&new_array[i].posX, sizeof(double), 1, fp2);
      fwrite(&new_array[i].posY, sizeof(double), 1, fp2);
      fwrite(&new_array[i].mass, sizeof(double), 1, fp2);
      fwrite(&new_array[i].velX, sizeof(double), 1, fp2);
      fwrite(&new_array[i].velY, sizeof(double), 1, fp2);
      fwrite(&new_array[i].b, sizeof(double), 1, fp2);
    }


    fclose(fp2);
    free(array);
    free(new_array);

  return 0;
}