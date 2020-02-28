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
  double centerOfMassX;
  double centerOfMassY;
} particleBox;

particle * read_particle(int N, FILE *fp1) {
  unsigned char buffer[8];
  double arr[6];
  particle *array = malloc(N * sizeof(particle));
  for(int i = 0; i<N; i++){
    for(int j = 0; j<6; j++){
        int e = fread(buffer, sizeof(buffer), 1, fp1);
        if(!e){
          printf("Error reading file!\n");
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
  if (box != NULL) printf("Mass: %f \n", (*box).mass);
  if ((*box).ne != NULL)
    print_stars((*box).ne);
  if ((*box).nw != NULL)
    print_stars((*box).nw);
  if ((*box).se != NULL)
    print_stars((*box).se);
  if ((*box).sw != NULL)
    print_stars((*box).sw);
}

void fitParticle(particleBox** box, particle star){
  if(*box == NULL){
    printf("Adding node \n");
    *box = (particleBox*)malloc(sizeof(particleBox));
    (**box).mass = star.mass;
    (**box).star = &star;
  }
  else if((**box).star == NULL){
    printf("Found parent node\n");
    (**box).mass += star.mass;
    //(**box).centerOfMassX = abs((**box).centerOfMassX-star.posX)*0.5;
    //(**box).centerOfMassY = abs((**box).centerOfMassX-star.posX)*0.5;
    printf("Star: x: %f y: %f\n", star.posX, star.posY);
    particleBox tempBox;
    if(star.posX <= (**box).x && star.posY > (**box).y){
      printf("nw\n");
      fitParticle(&(**box).nw, star);    
      (**box).nw = &tempBox;
      (*box)->nw->x = (**box).x*0.5;
      (*box)->nw->y = (**box).y*1.5;

      }else if(star.posX > (**box).x && star.posY > (**box).y){
        printf("ne\n");
        printf("x: %f y: %f \n", (**box).x, (**box).y);
        fitParticle(&(**box).ne, star);
        (**box).ne = &tempBox;
        (*box)->ne->x = (**box).x*1.5;
        (*box)->ne->y = (**box).y*1.5;

      }else if(star.posX <= (**box).x && star.posY <= (**box).y){
        printf("sw\n");
        fitParticle(&(**box).sw, star);
        (*box)->se->x = (**box).x*0.5;
        (*box)->sw->y = (**box).y*0.5;
      }else if(star.posX > (**box).x && star.posY <= (**box).y){
        printf("se\n");
        fitParticle(&(**box).se, star);
        (*box)->se->x = (**box).x*1.5;
        (*box)->se->y = (**box).y*0.5;
      }
      printf("Box: x: %f y: %f \n", tempBox.x, tempBox.y);
  }else{
    particle *oldStar = (**box).star;
    (**box).star = NULL;
    fitParticle(box, star);
    fitParticle(box, *oldStar);
    (**box).mass += star.mass;
  }
}

particleBox calcForce(particle star, particleBox box){
  if(box.star != NULL){
    
  }
}

// ./galsim 2 input_data/circles_N_2.gal 500 1e-5 0.1 0
int main(int argc, char* argv[]){  
  if (argc != 7){
      printf("Wrong number of input arguments\n");
      return 1;
  }  
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
  particle *array = read_particle(N, fp1);

  for(int i = 0; i<N; i++){
    printf("Particle %d: xpos: %f, ypos: %f, mass: %f \n", i, array[i].posX, array[i].posY, array[i].mass);
  }

  particleBox *root = (particleBox*)malloc(sizeof(particleBox));
  
  (*root).mass = 0;
  (*root).x = 0.5;
  (*root).y = 0.5;

  (*root).centerOfMassX = 0.5;
  (*root).centerOfMassY = 0.5;

  for(int i =0; i<2; i++){
    fitParticle(&root, array[i]);
  }
  print_stars(root);
  return 0;
}