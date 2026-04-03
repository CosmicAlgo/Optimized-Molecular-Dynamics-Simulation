/*
 * Header file for MD simulation
 * Optimized version for Coursework 2
 * 
 * This file defines constants and declarations for the
 * molecular dynamics simulation.
 */

#ifdef DECL
#define DEF
#else
#define DEF extern
#endif

/* 
 * System parameters
 * Nbody: Number of particles (4096 = 4 * 1024)
 * Npair: Number of unique particle pairs for pairwise interactions
 */
#define Nbody (4 * 1024)
#define Npair ((Nbody * (Nbody - 1)) / 2)

/* Dimension enumeration for array indexing */
enum { Xcoord = 0, Ycoord = 1, Zcoord = 2, Ndim = 3 };

/* 
 * Global arrays - allocated dynamically in control.c
 * Using contiguous memory layout for cache efficiency
 */
DEF double *pos[Ndim];      /* Particle positions [dim][particle] */
DEF double *velo[Ndim];     /* Particle velocities [dim][particle] */
DEF double *f[Ndim];        /* Forces on particles [dim][particle] */
DEF double *vis;           /* Viscosity coefficients per particle */
DEF double *mass;           /* Mass per particle */
DEF double *radius;         /* Radius per particle (for collision) */
DEF double *delta_pos[Ndim]; /* Pairwise separation vectors [dim][pair] */
DEF double *r;              /* Distance from central mass per particle */
DEF double *delta_r;        /* Pairwise separation distances */
DEF double wind[Ndim];      /* Wind vector (external force field) */
DEF int collisions;         /* Collision counter */

/* Physical constants */
#define G 2.0               /* Gravitational constant */
#define M_central 1000.0    /* Central mass for orbital dynamics */

/* Function declarations */
void evolve(int Nstep, double dt);

/* Trajectory output for visualization */
void write_trajectory_header(int sample_interval);
void write_trajectory_frame(void);
void close_trajectory(void);
