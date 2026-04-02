/*
 *  Utility functions for MD simulation.
 *  Optimized version for Coursework 2
 */
#include <math.h>

/*
 * Viscosity forces: f[i] = -vis[i] * velo[i]
 * Simple element-wise operation - compiler auto-vectorizes well
 */
void vis_forces(int N, double *f, double *vis, double *velo) {
  int i;
  for(i = 0; i < N; i++) {
    f[i] = -vis[i] * velo[i];
  }
}

/*
 * Wind forces: f[i] = f[i] - vis[i] * wind
 * Accumulates wind effect onto existing forces
 */
void wind_forces(int N, double *f, double *vis, double wind) {
  int i;
  for(i = 0; i < N; i++) {
    f[i] -= vis[i] * wind;
  }
}

/*
 * Add squared components: r[k] += delta[k]^2
 * Used for computing vector norms
 */
void add_norms(int N, double *r, double *delta) {
  int k;
  for(k = 0; k < N; k++) {
    r[k] += delta[k] * delta[k];
  }
}

/*
 * Force calculation: F = W * delta / r^3
 * 
 * NOTE: This function is retained for standalone use.
 * The main simulation in MD.c uses an inlined version
 * for performance in the hot loop.
 *
 * Original used pow(rv, 3.0) which is a general-purpose
 * library call. Direct multiplication rv*rv*rv is used instead.
 */
double forces(double Wv, double deltav, double rv) {
  double r3 = rv * rv * rv;
  return Wv * deltav / r3;
}
