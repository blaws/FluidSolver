// Fluid Solver
// blaws, 6/7/2013
// Based on "Real-Time Fluid Dynamics for Games", by Jos Stam

#define N 150
#define SIZE (N+2)*(N+2)
#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}

#include <GLUT/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static float* u;
static float* v;
static float* u_prev;
static float* v_prev;
static float* dens[3];
static float* dens_prev[3];
float visc = 10;
float dt = .1;
float diff = .1;
float xmult=3,ymult=3;
int currentButton=0,currentColor=1;

void add_source(float* x,float* s,float dt){
  int i;
  for(i=0;i<SIZE;i++) x[i]+=dt*s[i];
}

void set_bnd(int b,float* x){
  int i,m;

  for(m=0;m<3;m++){
    for(i=1;i<=N;i++){
      x[IX(0,i)] = (b==1) ? -x[IX(1,i)] : x[IX(1,i)];
      x[IX(N+1,i)] = (b==1) ? -x[IX(N,i)] : x[IX(N,i)];
      x[IX(i,0)] = (b==2) ? -x[IX(i,1)] : x[IX(i,1)];
      x[IX(i,N+1)] = (b==2) ? -x[IX(i,N)] : x[IX(i,N)];
    }

    x[IX(0,0)] = .5*(x[IX(1,0)]+x[IX(0,1)]);
    x[IX(0,N+1)] = .5*(x[IX(1,N+1)]+x[IX(0,N)]);
    x[IX(N+1,0)] = .5*(x[IX(N,0)]+x[IX(N+1,1)]);
    x[IX(N+1,N+1)] = .5*(x[IX(N,N+1)]+x[IX(N+1,N)]);
  }
}

void diffuse(int b,float* x,float* x0,float diff,float dt){
  int i,j,k;
  float a = dt * diff * N*N;

  for(k=0;k<20;k++){
    for(i=1;i<=N;i++){
      for(j=1;j<=N;j++){
	x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/(1+4*a);
      }
    }
    set_bnd(b,x);
  }
}

void advect(int b,float* d,float* d0,float* u,float* v,float dt){
  int i,j,i0,j0,i1,j1;
  float x,y,s0,t0,s1,t1,dt0;

  dt0 = dt*N;
  for(i=1;i<=N;i++){
    for(j=1;j<=N;j++){
      x = i-dt0*u[IX(i,j)];
      y = j-dt0*v[IX(i,j)];
      if(x<.5) x=.5;
      else if(x>N+.5) x=N+.5;
      i0 = (int)x;
      i1 = i0 + 1;
      if(y<.5) y=.5;
      else if(y>N+.5) y=N+.5;
      j0 = (int)y;
      j1 = j0 + 1;
      s1 = x - i0;
      s0 = 1 - s1;
      t1 = y - j0;
      t0 = 1 - t1;
      d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)]) + s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
    }
  }
  set_bnd(b,d);
}

void dens_step(float* x[],float* x0[],float* u,float* v,float diff,float dt){
  int m;
  for(m=0;m<3;m++){
    add_source(x[m],x0[m],dt);
    SWAP(x0[m],x[m]);
    diffuse(0,x[m],x0[m],diff,dt);
    SWAP(x0[m],x[m]);
    advect(0,x[m],x0[m],u,v,dt);
  }
}

void project(float* u,float* v,float* p,float* div){
  int i,j,k;
  float h;

  h = 1.0/N;
  for(i=1;i<=N;i++){
    for(j=1;j<=N;j++){
      div[IX(i,j)] = -.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]);
      p[IX(i,j)] = 0;
    }
  }
  set_bnd(0,div);
  set_bnd(0,p);

  for(k=0;k<20;k++){
    for(i=1;i<=N;i++){
      for(j=1;j<=N;j++){
	p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+p[IX(i,j-1)]+p[IX(i,j+1)])/4;
      }
    }
    set_bnd(0,p);
  }

  for(i=1;i<=N;i++){
    for(j=1;j<=N;j++){
      u[IX(i,j)] -= .5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
      v[IX(i,j)] -= .5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
    }
  }
  set_bnd(1,u);
  set_bnd(2,v);
}

void vel_step(float* u,float* v,float* u0,float* v0,float visc,float dt){
  add_source(u,u0,dt);
  add_source(v,v0,dt);
  SWAP(u0,u);
  diffuse(1,u,u0,visc,dt);
  SWAP(v0,v);
  diffuse(2,v,v0,visc,dt);
  project(u,v,u0,v0);
  SWAP(u0,u);
  SWAP(v0,v);
  advect(1,u,u0,u0,v0,dt);
  advect(2,v,v0,u0,v0,dt);
  project(u,v,u0,v0);
}


void draw_dens(void){
  //printf("Drawing... ");
  int x,y,m;//,width;
  float c[3];
  glClear(GL_COLOR_BUFFER_BIT);

  for(y=0;y<=N;y++){
    for(x=0;x<=N;x++){
      //width = 0;
      //while(x+width<=N && dens[m][IX(x,y)]==dens[m][IX(x+width,y)]) width++;
      for(m=0;m<3;m++) c[m] = dens[m][IX(x,y)]/255.0;
      if(c[0]==0 && c[1]==0 && c[2]==0) continue;
      glColor3f(c[0],c[1],c[2]);

      glRecti(x,N-y,x+1,N-y+1);
      //x+=width-1;
    }
  }
  glutSwapBuffers();
  //printf("done\n");
}

void reshape(int w,int h){
  glViewport(0,0,(GLsizei)w,(GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0,N,0,N,-1,1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  xmult = w/(GLfloat)N;
  ymult = h/(GLfloat)N;
  //printf("xmult=%f\nymult=%f\n",xmult,ymult);
}

void mouse(int button,int state,int x,int y){
  int i,j,m,newx=x/xmult,newy=y/ymult;
  switch(button){
  case GLUT_LEFT_BUTTON:
    if(state==GLUT_DOWN){
      currentButton=1;
      for(i=(newx-5>=0?newx-5:0); i<newx+5 && i<=N; i++){
	for(j=(newy-5>=0?newy-5:0); j<newy+5 && j<=N; j++){
	  dens_prev[currentColor][IX(i,j)] = 255.0;
	}
      }
      glutPostRedisplay();
      //printf("Left mouse down (%d,%d).\n",newx,newy);
    }
    else{
      currentButton=0;
      currentColor = (currentColor+1)%3;
    }
    break;
  case GLUT_RIGHT_BUTTON:
    if(state==GLUT_DOWN){
      currentButton=2;
      for(i=(newx-5>=0?newx-5:0); i<newx+5 && i<=N; i++){
	for(j=(newy-5>=0?newy-5:0); j<newy+5 && j<=N; j++){
	  for(m=0;m<3;m++) dens_prev[m][IX(i,j)] = 0.0;
	}
      }
      glutPostRedisplay();
      //printf("Right mouse down (%d,%d).\n",newx,newy);
    }
    else currentButton=0;
    break;
  case GLUT_MIDDLE_BUTTON:
    if(state==GLUT_DOWN){
      currentButton=3;
      for(i=0;i<=N;i++){
	for(j=0;j<=N;j++){
	  //u_prev[IX(i,j)] += 10.0;
	  v_prev[IX(i,j)] += 10.0;
	}
      }
      glutPostRedisplay();
      //printf("Middle mouse down (%d,%d).\n",newx,newy);
    }
    else currentButton=0;
    break;
  default:
    break;
  }
}

void mouseMove(int x,int y){
  int i,j,m,newx=x/xmult,newy=y/ymult;
  if(currentButton==3){
    for(i=0;i<=N;i++){
      for(j=0;j<=N;j++){
	v_prev[IX(i,j)] = 10;
      }
    }
  }
  else{
    for(i=(newx-5>=0?newx-5:0); i<newx+5 && i<=N; i++){
      for(j=(newy-5>=0?newy-5:0); j<newy+5 && j<=N; j++){
	if(currentButton==1) dens_prev[currentColor][IX(i,j)] = 255.0;
	else if(currentButton==2) for(m=0;m<3;m++) dens_prev[m][IX(i,j)] = 0.0;
      }
    }
  }
  glutPostRedisplay();
}

void fluidMainLoop(void){
  //printf("Loop... ");
  vel_step(u,v,u_prev,v_prev,visc,dt);
  dens_step(dens,dens_prev,u,v,diff,dt);
  glutPostRedisplay();
  //printf("done.\n");
}

int main(int argc,char* argv[]){
  int i,j;
  u = malloc(SIZE*sizeof(float));
  v = malloc(SIZE*sizeof(float));
  u_prev = malloc(SIZE*sizeof(float));
  v_prev = malloc(SIZE*sizeof(float));
  for(j=0;j<3;j++){
    dens[j] = malloc(SIZE*sizeof(float));
    dens_prev[j] = malloc(SIZE*sizeof(float));
  }
  if(dens[0]==NULL || dens[1]==NULL || dens[2]==NULL || dens_prev[0]==NULL || dens_prev[1]==NULL || dens_prev[2]==NULL){
    printf("Memory allocation failure.\n");
    exit(1);
  }

  for(i=0;i<SIZE;i++){
    u[i] = 0.0;
    v[i] = 0.0;
    u_prev[i] = 0.0;
    v_prev[i] = 0.0;
    for(j=0;j<3;j++){
      dens[j][i] = 0.0;
      dens_prev[j][i] = 0.0;
    }
  }

  // GLUT
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
  glutInitWindowSize(3*N,3*N);
  glutInitWindowPosition(0,0);
  glutCreateWindow(argv[0]);
  glClearColor(0,0,0,0);
  glShadeModel(GL_FLAT);
  glOrtho(0,N,0,N,-1,1);
  glutDisplayFunc(draw_dens);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMove);
  glutIdleFunc(fluidMainLoop);
  glutMainLoop();


  free(u);
  free(v);
  free(u_prev);
  free(v_prev);
  free(dens);
  free(dens_prev);
}
