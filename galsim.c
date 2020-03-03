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
  double Fx, Fy;
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

    //box->mass += star->mass;
    //printf("Mass for box %f, %f: %f\n",box->x, box->y, box->mass);

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
    //box->mass += star->mass;
    //printf("Mass for box %f, %f: %f\n",box->x, box->y, box->mass);

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


void calcCenterOfMass(particleBox * box){
  double r1, r2, r3, r4;
  double m1, m2, m3, m4;
  double xCenter, yCenter;

  r1 = box->nw->centerOfMassX;
  r2 = box->ne->centerOfMassX;
  r3 = box->sw->centerOfMassX;
  r4 = box->se->centerOfMassX;

  m1 = box->nw->mass;
  m2 = box->ne->mass;
  m3 = box->sw->mass;
  m4 = box->se->mass;

  xCenter = (r1*m2 + r2*m2 + r3*m3 + r4*m4)/(m1 + m2 + m3 + m4);

  r1 = box->nw->centerOfMassY;
  r2 = box->ne->centerOfMassY;
  r3 = box->sw->centerOfMassY;
  r4 = box->se->centerOfMassY;

  yCenter = (r1*m2 + r2*m2 + r3*m3 + r4*m4)/(m1 + m2 + m3 + m4);

  box->centerOfMassX = xCenter;
  box->centerOfMassY = yCenter;
}

void calcMass(particleBox * box){
  if(box->isEmpty){
    box->mass = 0;
    box->centerOfMassX = 0;
    box->centerOfMassY = 0;

  }else if(box->isLeaf){
    box->mass=box->star->mass;
    box->centerOfMassX = box->star->posX;
    box->centerOfMassY = box->star->posY;

  }else{
    calcMass(box->nw);
    calcMass(box->ne);
    calcMass(box->sw);
    calcMass(box->se);

    box->mass += box->nw->mass;
    box->mass += box->ne->mass;
    box->mass += box->sw->mass;
    box->mass += box->se->mass;

    calcCenterOfMass(box);
  }
}

void calcForce(particle * star, particleBox * box, double thetaMax){
  double rx, ry;
  double e0 = 0.001;
  rx = abs(star->posX - box->centerOfMassX);
  ry = abs(star->posY - box->centerOfMassY);
  double r = sqrt(rx*rx+ry*ry);
  if(box->star != star){
    if(box->isEmpty){
      
    }
    else if(box->isLeaf){
      double denom = (r+e0)*(r+e0)*(r+e0);
      star->Fx += star->mass * rx / denom;
      star->Fy += star->mass * ry / denom;
    } 
    else if (box->side / r <= thetaMax) {
      double denom = (r+e0)*(r+e0)*(r+e0);
      star->Fx += star->mass * rx / denom;
      star->Fy += star->mass * ry / denom;
    } else {
      calcForce(star, box->nw, thetaMax);
      calcForce(star, box->ne, thetaMax);
      calcForce(star, box->sw, thetaMax);
      calcForce(star, box->se, thetaMax);
    }
  }
}

void deleteBoxes(particleBox * box) {
  if (box->isLeaf) {
    free(box);
    box = NULL;
  } else {
    deleteBoxes(box->nw);
    deleteBoxes(box->ne);
    deleteBoxes(box->sw);
    deleteBoxes(box->se);
    free(box);
    box = NULL;
  }
}

// ./galsim 2 input_data/circles_N_2.gal 500 1e-5 0.1 0
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

  printf("Command line arguments given: %d, %s, %d, %f, %f, %d \n", N, filename, n_steps, delta_t, theta_max, graphics);
  const double G = 100.0 / N;
  


  particle **array = read_particle(N, filename);

  for(int i = 0; i<n_steps; i++){
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

    calcMass(root);
    for(int i =0; i<N; i++){
      array[i]->Fx = 0;
      array[i]->Fy = 0;

      calcForce(array[i], root, theta_max);

      array[i]->Fx *= -G*array[i]->mass;
      array[i]->Fy *= -G*array[i]->mass;

    }
    deleteBoxes(root);

    //print_stars(root);
  }
  return 0;
}

//./galsim 4 input_data/circles_N_4.gal 500 1e-5 0.1 0