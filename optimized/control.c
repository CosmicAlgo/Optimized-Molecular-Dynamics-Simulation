/*
 *
 * Control program for the MD update
 * Optimized version for Coursework 2
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define DECL
#include "coord.h"

double second(void);

int main(int argc, char *argv[]) {
  int i, j;
  FILE *in, *out;
  double tstart, tstop;
  double start, stop;
  char name[80];
  
  /* Timestep value - fixed for simulation stability */
  double dt = 0.02;
  
  /* Number of timesteps to use */
  int Nstep = 100;
  int Nsave = 5;
  
  if(argc > 1) {
    Nstep = atoi(argv[1]);
  }
  
  /* Wind vector - external force field */
  wind[Xcoord] = 0.9;
  wind[Ycoord] = 0.4;
  wind[Zcoord] = 0.0;
  
  /* 
   * Memory allocation for 3D arrays
   * Using contiguous allocation for better cache locality
   */
  r = calloc(Nbody, sizeof(double));
  mass = calloc(Nbody, sizeof(double));
  radius = calloc(Nbody, sizeof(double));
  vis = calloc(Nbody, sizeof(double));
  
  /* Contiguous allocation for 2D arrays.
   * Single calloc per array type with pointer arithmetic for rows.
   * This guarantees all dimensions are contiguous in memory,
   * improving cache prefetch behaviour.
   *
   * Note: delta_pos[] and delta_r[] are no longer allocated.
   * The optimised evolve() computes pairwise separations on-the-fly
   * using register variables, eliminating ~192MB of temporary storage
   * (Ndim * Npair * sizeof(double) = 3 * 8386560 * 8 bytes). */
  f[0] = calloc(Ndim * Nbody, sizeof(double));
  pos[0] = calloc(Ndim * Nbody, sizeof(double));
  velo[0] = calloc(Ndim * Nbody, sizeof(double));
  
  for(i = 1; i < Ndim; i++) {
    f[i] = f[0] + i * Nbody;
    pos[i] = pos[0] + i * Nbody;
    velo[i] = velo[0] + i * Nbody;
  }

  /* Initialize collision counter */
  collisions = 0;
  
  /* Read initial particle data from file */
  in = fopen("input.dat", "r");
  if(!in) {
    perror("input.dat");
    exit(1);
  }

  /* Read: mass, radius, viscosity, position (x,y,z), velocity (x,y,z) */
  for(i = 0; i < Nbody; i++) {
    fscanf(in, "%16le%16le%16le%16le%16le%16le%16le%16le%16le\n",
      mass + i, radius + i, vis + i,
      &pos[Xcoord][i], &pos[Ycoord][i], &pos[Zcoord][i],
      &velo[Xcoord][i], &velo[Ycoord][i], &velo[Zcoord][i]);
  }
  fclose(in);

  /* 
   * Run simulation and measure performance
   * Timing excludes file I/O as per specification
   */
  tstart = second();
  
  for(j = 1; j <= Nsave; j++) {
    start = second();
    evolve(Nstep, dt);
    stop = second();
    
    printf("%d timesteps took %f seconds\n", Nstep, stop - start);
    printf("collisions %d\n", collisions);
    fflush(stdout);
    
    /* Write checkpoint output */
    sprintf(name, "output.dat%03d", j * Nstep);
    out = fopen(name, "w");
    if(!out) {
      perror(name);
      exit(1);
    }

    for(i = 0; i < Nbody; i++) {
      fprintf(out, "%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E\n",
        mass[i], radius[i], vis[i],
        pos[Xcoord][i], pos[Ycoord][i], pos[Zcoord][i],
        velo[Xcoord][i], velo[Ycoord][i], velo[Zcoord][i]);
    }
    fclose(out);
  }
  
  tstop = second();
  printf("Total: %d timesteps took %f seconds\n", Nsave * Nstep, tstop - tstart);

  return 0;
}

/* 
 * High-resolution timer using gettimeofday
 * Returns time in seconds with microsecond precision
 */
double second() {
  struct timeval tp;
  struct timezone tzp;
  
  gettimeofday(&tp, &tzp);
  return (double)tp.tv_sec + (double)tp.tv_usec * 1.e-6;
}
