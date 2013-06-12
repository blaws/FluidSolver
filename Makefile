all: FluidSolver

FluidSolver: FluidSolver.c
	gcc FluidSolver.c -o FluidSolver -framework GLUT -framework OpenGL

clean:
	rm -f FluidSolver *~
