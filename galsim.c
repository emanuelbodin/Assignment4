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
  double Fx, Fy;
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

void writeToFile(particle ** array, int N) {
  char * filename = "result.gal";
  FILE * fp;
  fp = fopen(filename, "wb");
  for (int i = 0; i < N; i++) {
    fwrite(array[i], sizeof(double)*6, 1, fp);
  }
  fclose(fp);
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
    else {
      printf("Error when trying to insert particle");
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

  xCenter = (r1*m1 + r2*m2 + r3*m3 + r4*m4)/(m1 + m2 + m3 + m4);

  r1 = box->nw->centerOfMassY;
  r2 = box->ne->centerOfMassY;
  r3 = box->sw->centerOfMassY;
  r4 = box->se->centerOfMassY;

  yCenter = (r1*m1 + r2*m2 + r3*m3 + r4*m4)/(m1 + m2 + m3 + m4);

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
  double rx, ry, denom;
  double e0 = 0.001;
  rx = star->posX - box->centerOfMassX;
  ry = star->posY - box->centerOfMassY;
  double r = sqrt(rx*rx+ry*ry);
 // if(box->star != star){
/*     if(box->isEmpty){
      
    } */
     if(box->isLeaf){
      denom = (r+e0)*(r+e0)*(r+e0);
      star->Fx += star->mass * rx / denom;
      star->Fy += star->mass * ry / denom;
    } 
    else if ((box->side / r) < thetaMax) {
      denom = (r+e0)*(r+e0)*(r+e0);
      star->Fx += star->mass * rx / denom;
      star->Fy += star->mass * ry / denom;
    } else if ((box->side / r) > thetaMax) {
      calcForce(star, box->nw, thetaMax);
      calcForce(star, box->ne, thetaMax);
      calcForce(star, box->sw, thetaMax);
      calcForce(star, box->se, thetaMax);
    } else {
      printf("Error when calculating force \n");
    }
  //}
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

void printArray(particle ** a, int N) {
  for (int i = 0; i < N; i++) {
    printf("mass: %f, Xpos: %f, Ypos: %f, b: %f \n", a[i]->mass, a[i]->posX, a[i]->posY, a[i]->b);
  }
}

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
  printArray(array, N);
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
    for(int j = 0; j<N; j++){
      fitParticle(root, array[j]);
    }

    calcMass(root);
    for(int j =0; j<N; j++){
      array[j]->Fx = 0;
      array[j]->Fy = 0;

      calcForce(array[j], root, theta_max);

      array[j]->Fx *= -G;
      array[j]->Fy *= -G;

    }
    deleteBoxes(root);

    for (int j = 0; j < N; j++) {
      array[j]->velX += delta_t*(array[j]->Fx);
      array[j]->velY += delta_t*(array[j]->Fy);
      array[j]->posX += delta_t*array[j]->velX;
      array[j]->posY += delta_t*array[j]->velY;
    }
  }
  printf("\n");
  printArray(array, N);
  writeToFile(array, N);
  for (int i = 0; i < N; i++) {
    free(array[i]);
    array[i] = NULL;
  }
  free(array);
  array = NULL;
  return 0;
}

//./galsim 4 input_data/circles_N_4.gal 500 1e-5 0.1 0
//./galsim 10 input_data/ellipse_N_00010.gal 200 1e-5 0.1 0

//./compare_gal_files 10 ../result.gal ../ref_output_data/ellipse_N_00010_after200steps.gal