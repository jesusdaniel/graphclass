/*********************************************************************
  **********************************************************************
  ** TO COMPILE:
  ** gcc -shared -fPIC -o Pmax2.so Pmax2.c -lm
** R CMD SHLIB foo.c    OR    R --arch x64 CMD SHLIB foo.c
*********************************************************************/
  # include <stdlib.h>
  # include <stdio.h>
  # include <math.h>
  # include <limits.h>
  # include <float.h>
  
  # define EPS DBL_EPSILON
  
  void pmax2(double *x, int *d){
    int i;
    for (i=0; i<d[0]; i++)
      if(x[i] < 0)
        x[i] = 0;
  };

void soft_thresholding(double *x, double *lambda, int *d) {
  int i;
  for(i=0; i<d[0]; i++){
    if(x[i]>lambda[0]){
      x[i] = x[i]-lambda[0];
    }else{
      if(x[i]< -lambda[0]){
        x[i] = x[i] + lambda[0];
      }
      else{
        x[i] = 0;
      }
    }
  }
};

void group_lasso_node_penalty(double *Db, int *N, double *gl){
  int i,j;
  double norm_node;
  gl[0] = 0;
  for(i=0; i <N[0]; i++) {
    norm_node = 0;
    for(j=(N[0]-1)*i; j<(N[0]-1)*(i+1); j++) {
      norm_node += Db[j]*Db[j];
    }
    gl[0] += sqrt(norm_node);
  }
}
void l1l2_node_soft_thresholding(double *b, int *N, double *lambda1, double * lambda2) {
  int i,j,p;
  p = N[0]*(N[0]-1);
  for(i=0; i<p; i++){
    if(b[i]>lambda1[0]){
      b[i] = b[i]-lambda1[0];
    }else{
      if(b[i]< -lambda1[0]){
        b[i] = b[i] + lambda1[0];
      }
      else{
        b[i] = 0;
      }
    }
  }
  double norm_node, t;
  for(i=0; i <N[0]; i++) {
    norm_node = 0;
    for(j=(N[0]-1)*i; j<(N[0]-1)*(i+1); j++) {
      norm_node += b[j]*b[j];
    }
    norm_node = sqrt(norm_node);
    if(norm_node<=lambda2[0]) {
      for(j=(N[0]-1)*i; j<(N[0]-1)*(i+1); j++){
        b[j] = 0;
      }
    }
    else{
      t = 1-lambda2[0]/norm_node;
      for(j=(N[0]-1)*i; j<(N[0]-1)*(i+1); j++){
        b[j] = b[j]*t;
      }
    }
  }
}
void l2_node_soft_ghtesholding(double *Db, int *N, double *lambda){
  int i,j;
  double norm_node;
  double t;
  for(i=0; i <N[0]; i++) {
    norm_node = 0;
    for(j=(N[0]-1)*i; j<(N[0]-1)*(i+1); j++) {
      norm_node += Db[j]*Db[j];
    }
    norm_node = sqrt(norm_node);
    if(norm_node<=lambda[0]) {
      for(j=(N[0]-1)*i; j<(N[0]-1)*(i+1); j++){
        Db[j] = 0;
      }
    }
    else{
      t = 1-lambda[0]/norm_node;
      for(j=(N[0]-1)*i; j<(N[0]-1)*(i+1); j++){
        Db[j] = Db[j]*t;
      }
    }
  }
}

void l2_soft_thresholding(double *x, double *lambda, int *d) {
  int i;
  double norm;
  norm = 0;
  for(i=0; i<d[0]; i++){
    norm += x[i]*x[i];
  }
  norm = sqrt(norm);
  if(norm <= lambda[0]) {
    for(i=0; i<d[0]; i++){
      x[i] = 0;
    }
  }
  else{
    double t;
    t = 1-lambda[0]/norm;
    for(i=0; i<d[0]; i++){
      x[i] = x[i]*t;
    }
  }
};
