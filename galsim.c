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
  double side;

  double mass;
  double centerOfMassX;
  double centerOfMassY;

  int isEmpty;
  int isLeaf;
} particleBox;

particle ** read_particle(int N, char * filename) {
  FILE* file = fopen(filename, "rb");
  double * buffer = (double*)malloc(sizeof(double)*6*N);

  int readOk = fread(buffer, sizeof(double)*6*N, 1, file);

  fclose(file);

  if(readOk != 1){
    printf("Error reading file");
    return NULL;
  }
  particle **array = (particle**)malloc(N * sizeof(particle*));


  for(int i = 0; i<N; i++){
    array[i] = (particle*)malloc(sizeof(particle));
  }
  int offset = 0;
  for(int i = 0; i<N; i++){
    array[i]->posX = *(buffer+offset);
    offset++;
    array[i]->posY = *(buffer+offset);
    offset++;
    array[i]->mass = *(buffer+offset);
    offset++;
    array[i]->velX = *(buffer+offset);
    offset++;
    array[i]->velY = *(buffer+offset);
    offset++;
    array[i]->b = *(buffer+offset);
    offset++;
  }

  free(buffer);
  buffer = NULL;
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

particleBox * createBox(double x, double y, double side){
  particleBox *box;

  box = malloc(sizeof(particleBox));

  box->mass = 0;
  box->centerOfMassX = 0;
  box->centerOfMassY = 0;

  box->x = x;
  box->y = y;
  box->side = side;

  box->star = NULL;
  box->isEmpty = 1;
  box->isLeaf = 1;


  return box;
}

void fitParticle(particleBox * box, particle * star){
  if(box->isEmpty == 1){
    box->star = star;
    box->isEmpty = 0;
    box->isLeaf = 1;
  }
  else if(box->isLeaf == 1){
    particle *oldStar = box->star;
    box->star = NULL;
    box->isLeaf = 0;

    box->nw = createBox(box->x-box->side/4, box->y+box->side/4, box->side/2);
    box->ne = createBox(box->x+box->side/4, box->y+box->side/4, box->side/2);
    box->sw = createBox(box->x-box->side/4, box->y-box->side/4, box->side/2);    
    box->se = createBox(box->x+box->side/4, box->y-box->side/4, box->side/2);

    fitParticle(box, star);
    fitParticle(box, oldStar);
    //box->mass += star.mass;
      
  }else{

    //box->mass += star.mass;
    //box->centerOfMassX = abs(box->centerOfMassX-star.posX)*0.5;
    //box->centerOfMassY = abs(box->centerOfMassX-star.posX)*0.5;
    if(star->posX <= box->x && star->posY > box->y){
      fitParticle(box->nw, star);
      }else if(star->posX > box->x && star->posY > box->y){
        fitParticle(box->ne, star);
      }else if(star->posX <= box->x && star->posY <= box->y){
        fitParticle(box->sw, star);
      }else if(star->posX > box->x && star->posY <= box->y){
        fitParticle(box->se, star);
      }
  }
}

//particleBox calcForce(particle star, particleBox box){
//  if(box.star != NULL){
//    
//  }
//}

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

 int N = atoi(argv[1]);
 char* filename = argv[2];
  particle **array = read_particle(N, filename);


  particleBox *root = (particleBox*)malloc(sizeof(particleBox));
  
  (*root).mass = 0;
  (*root).x = 0.5;
  (*root).y = 0.5;
  (*root).side = 1;
  (*root).star = NULL;

  (*root).nw = NULL;
  (*root).nw = NULL;
  (*root).nw = NULL;
  (*root).nw = NULL;

  (*root).isEmpty = 1;
  (*root).isLeaf = 1;

  (*root).centerOfMassX = 0.5;
  (*root).centerOfMassY = 0.5;

  for(int i =0; i<N; i++){
    //printf("Particle %d: xpos: %f, ypos: %f, mass: %f \n", i, array[i]->posX, array[i]->posY, array[i]->mass);
    fitParticle(root, array[i]);
  }
  //print_stars(root);
  return 0;
}

//./galsim 4 input_data/circles_N_4.gal 500 1e-5 0.1 0