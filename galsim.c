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
  double ax, ay;
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
    

    box->nw = createBox(box->x - box->side/4, box->y + box->side/4, box->side/2);
    box->ne = createBox(box->x + box->side/4, box->y + box->side/4, box->side/2);
    box->sw = createBox(box->x - box->side/4, box->y - box->side/4, box->side/2);    
    box->se = createBox(box->x + box->side/4, box->y - box->side/4, box->side/2);

    box->star = NULL;
    box->isLeaf = 0;

    fitParticle(box, oldStar);

    fitParticle(box, star);
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

void insert_mass(particleBox *n) {
	if (n->isEmpty) {
		// do nothing
	} else if (n->isEmpty) {
		n->mass = n->star->mass;
		n->centerOfMassX = n->star->posX;
		n->centerOfMassY = n->star->posY;
	} else {
		insert_mass(n->nw);
		insert_mass(n->ne);
		insert_mass(n->sw);
		insert_mass(n->se);
		
		// The mass of each parent node is the sum of the masses of its children
		double mass = n->nw->mass + n->ne->mass + n->sw->mass + n->se->mass;
		n->mass += mass;
		
		// "x*m" for each child
		double pos_m_nw = n->nw->centerOfMassX*n->nw->mass;
		double pos_m_ne = n->ne->centerOfMassX*n->ne->mass;
		double pos_m_sw = n->sw->centerOfMassX*n->sw->mass;
		double pos_m_se = n->se->centerOfMassX*n->se->mass;
		n->centerOfMassX = (pos_m_nw + pos_m_ne + pos_m_sw + pos_m_se)/mass;
		
		pos_m_nw = n->nw->centerOfMassY*n->nw->mass;
		pos_m_ne = n->ne->centerOfMassY*n->ne->mass;
		pos_m_sw = n->sw->centerOfMassY*n->sw->mass;
		pos_m_se = n->se->centerOfMassY*n->se->mass;
		n->centerOfMassY = (pos_m_nw + pos_m_ne + pos_m_sw + pos_m_se)/mass;
		
	}
}

void calcMass(particleBox * box){
  if(box->isEmpty){
    //box->mass = 0;
    //box->centerOfMassX = 0;
    //box->centerOfMassY = 0;

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
  double e0 = 1e-3;
  rx = star->posX - box->centerOfMassX;
  ry = star->posY - box->centerOfMassY;
  double r = sqrt(rx*rx+ry*ry);

  //double dx = box->x-star->posX;
  //double dy = box->y-star->posY;
  //double d = sqrt(dx*dx+dy*dy);

  

     if(box->isLeaf || (box->side / r) <= thetaMax){
      denom = 1/((r+e0)*(r+e0)*(r+e0));
      star->ax += star->mass * rx*denom;
      star->ay += star->mass * ry*denom;
    } else if ((box->side / r) > thetaMax) {
      calcForce(star, box->nw, thetaMax);
      calcForce(star, box->ne, thetaMax);
      calcForce(star, box->sw, thetaMax);
      calcForce(star, box->se, thetaMax);
    } else {
      printf("Error when calculating force \n");
    }
}


void calc_force(particle *p, particleBox *n, double theta_max) {
	double rx, ry, r, denom;
	double e0 = 1e-3;
	rx = p->posX<-n->centerOfMassX;
	ry = p->posY-n->centerOfMassY;
	r = sqrt(rx*rx+ry*ry);
	if (n->isLeaf) {
		denom = 1/((r+e0)*(r+e0)*(r+e0));
		p->ax += n->mass*denom*rx;
		p->ay += n->mass*denom*ry;
	} else if ((n->side/r) < theta_max) {
		denom = 1/((r+e0)*(r+e0)*(r+e0));
		p->ax += n->mass*denom*rx;
		p->ay += n->mass*denom*ry;
	} else if ((n->side/r) > theta_max){
		calc_force(p, n->nw, theta_max);
		calc_force(p, n->ne, theta_max);
		calc_force(p, n->sw, theta_max);
		calc_force(p, n->se, theta_max);
	} else {
		printf("Error in calc_force\n");
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

void printArray(particle ** a, int N) {
  for (int i = 0; i < N; i++) {
    printf("mass: %f, Xpos: %f, Ypos: %f, b: %f \n", a[i]->mass, a[i]->posX, a[i]->posY, a[i]->b);
  }
}

void printTree(particleBox * box){
  if(box->isLeaf){
    printf("Leaf node: \n");
    printf("side: %f, posx: %f, posy: %f, COMx: %f, COMy: %f\n", box->side, box->x, box->y, box->centerOfMassX, box->centerOfMassY);
    if(box->isEmpty){
      printf("This box is empty\n");
    }else{
      printf("This node contains the star: posx: %f, posy: %f, mass: %f\n", box->star->posX, box->star->posY, box->star->mass);
    }
    printf("\n");
  }else{
    printf("Parent node\n");
    printf("side: %f, posx: %f, posy: %f, COMx: %f, COMy: %f\n", box->side, box->x, box->y, box->centerOfMassX, box->centerOfMassY);
    printf("Has children: \n");
    printf("\n");
    printTree(box->nw);
    printTree(box->ne);
    printTree(box->sw);
    printTree(box->se);
    printf("STOP\n");
  }
}

particleBox* makeTree(particle **array, int N){
  particleBox *root;
  root = malloc(sizeof(particleBox));
  
    root->mass = 0;
    root->x = 0.5;
    root->y = 0.5;
    root->side = 1;
    root->star = NULL;

    root->nw = NULL;
    root->ne = NULL;
    root->sw = NULL;
    root->se = NULL;

    root->isEmpty = 1;
    root->isLeaf = 1;

    root->centerOfMassX = 0;
    root->centerOfMassY = 0;

    

    for(int j = 0; j<N; j++){
      fitParticle(root, array[j]);
    }
    return root;
}

particle** runSimulation(particle **array, int n_steps, double delta_t, int N, double theta_max){
  double G = (double)100/N;
  particleBox *root = NULL;

  for(int i = 0; i<=n_steps; i++){
    root = makeTree(array, N);

    printf("\n");


    printf("timestep %d. BEFORE CALC MASS First star values: \n", i);
    printf("posX: %f, posY: %f, mass: %f velY: %f, velY: %f, ax: %f, ay: %f\n", array[1]->posX, array[1]->posY, array[1]->mass, array[1]->velX, array[1]->velY, array[1]->ax, array[1]->ay);
    printf("\n");

    calcMass(root);

    //printTree(root);
    printf("timestep %d. First star values: \n", i);
    printf("posX: %f, posY: %f, mass: %f velY: %f, velY: %f, ax: %f, ay: %f\n", array[1]->posX, array[1]->posY, array[1]->mass, array[1]->velX, array[1]->velY, array[1]->ax, array[1]->ay);
    printf("\n");
    for(int j =0; j<N; j++){
      array[j]->ax = 0;
      array[j]->ay = 0;

      calcForce(array[j], root, theta_max);
      printf("timestep %d. After calc force\n", i);

      printf("posX: %f, posY: %f, mass: %f velY: %f, velY: %f, ax: %f, ay: %f\n", array[1]->posX, array[1]->posY, array[1]->mass, array[1]->velX, array[1]->velY, array[1]->ax, array[1]->ay);

      printf("\n");

      array[j]->ax *= -G;
      array[j]->ay *= -G;

    }
    printf("timestep %d. After multiplication with G First star values: \n", i);
    printf("posX: %f, posY: %f, mass: %f velY: %f, velY: %f, ax: %f, ay: %f\n", array[1]->posX, array[1]->posY, array[1]->mass, array[1]->velX, array[1]->velY, array[1]->ax, array[1]->ay);
    printf("\n");
    deleteBoxes(root);

    for (int j = 0; j < N; j++) {
      array[j]->velX += delta_t*(array[j]->ax);
      array[j]->velY += delta_t*(array[j]->ay);
      array[j]->posX += delta_t*array[j]->velX;
      array[j]->posY += delta_t*array[j]->velY;
    }
  }
  return array;
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
  
  particle **array = read_particle(N, filename);
 
  printArray(array, N);

  array = runSimulation(array, n_steps, delta_t, N, theta_max);

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