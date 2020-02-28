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
} particle;

typedef struct _particleBox
{
  struct _particleBox *nw;
  struct _particleBox *ne;
  struct _particleBox *sw;
  struct _particleBox *se;

  particle *star;

  double x;
  double y;

  double mass;
} particleBox;

particle * read_particle(int N, FILE *fp1) {
  unsigned char buffer[8];
  double arr[6];
  particle *array = malloc(N * sizeof(particle));
  for(int i = 0; i<N; i++){
    for(int j = 0; j<6; j++){
        int e = fread(buffer, sizeof(buffer), 1, fp1);
        if(e){
          printf("Error reading file");
          return NULL;
        }
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

void print_stars(particleBox *box) {
  if (box == NULL) printf("No stars! \n");
  if (box != NULL) printf("%f \n", (*box).mass);
  if ((*box).ne != NULL)
    print_stars((*box).ne);
  if ((*box).nw != NULL)
    print_stars((*box).nw);
  if ((*box).se != NULL)
    print_stars((*box).se);
  if ((*box).sw != NULL)
    print_stars((*box).sw);
}

particleBox fitParticle(particleBox** box, particle star){

  if(box == NULL){
    particleBox newBox;
    newBox.mass = star.mass;
    newBox.star = &star;
    return newBox;
  }

  if((**box).star == NULL){
    (**box).mass += star.mass;

    if(star.posX < (**box).x && star.posY > (**box).y){
      particleBox tempBox = fitParticle(&(**box).nw, star);
      (**box).nw = &tempBox;
      (tempBox).x = (**box).x*0.5;
      (tempBox).y = (**box).y*1.5;

      }else if(star.posX > (**box).x && star.posY < (**box).y){
        particleBox tempBox = fitParticle(&(**box).ne, star);
        (**box).ne = &tempBox;
        (tempBox).x = (**box).x*1.5;
        (tempBox).y = (**box).y*1.5;

      }else if(star.posX < (**box).x && star.posY < (**box).y){
        particleBox tempBox = fitParticle(&(**box).sw, star);
        (**box).sw = &tempBox;
        (tempBox).x = (**box).x*0.5;
        (tempBox).y = (**box).y*0.5;

      }else if(star.posX > (**box).x && star.posY > (**box).y){
        particleBox tempBox = fitParticle(&(**box).se, star);
        (**box).se = &tempBox;
        (tempBox).x = (**box).x*1.5;
        (tempBox).y = (**box).y*0.5;
      }
  }else{
    particle *oldStar = (**box).star;
    (**box).star = NULL;
    fitParticle(box, star);
    fitParticle(box, *oldStar);
    (**box).mass += star.mass;
  }
  return **box;
}

// ./galsim 2 input/cicles_N_2.gal 500 1e-5 0.1 0
int main(int argc, char* argv[]){
  printf("adihfihfiw %d", argc);
  /*
  if (argc != 6){
      printf("Wrong number of input arguments\n");
      return 1;
  }*/
  printf("hej");
  /*
  
  
  int n_steps = atoi(argv[3]);
  double delta_t = atof(argv[4]);
  double theta_max = atof(argv[5]);
  int graphics = atoi(argv[6]);

  printf("Command line arguments given: %d, %s, %d, %f, %f, %d \n", N, filename, n_steps, delta_t, theta_max, graphics);
  const double e0 = 0.001;
  const double G = 100.0 / N;
  */
 FILE *fp1;
 int N = atoi(argv[1]);
 char* filename = argv[2];
  fp1 = fopen(filename, "rb");
  if (fp1 == NULL){
  printf("Error while opening the file.\n");
  return 1;
  }
  printf("%s", "hej");
  particle *array = read_particle(N, fp1);

  particleBox *root = NULL;
  (*root).mass = 0;
  (*root).x = 0.5;
  (*root).y = 0.5;
  printf("%s", "hej");
  for(int i =0; i<N; i++){
    fitParticle(&root, array[i]);
  }
  printf("%s", "hej");
  print_stars(root);

  return 0;
}