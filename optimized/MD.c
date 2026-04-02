/*
 *  Simple molecular dynamics code.
 *  Optimized version for Coursework 2 - Performance Programming
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"

/*
 * Note on forces() replacement:
 * The original util.c function  return Wv*deltav/pow(rv,3.0)  is not
 * called from this file. Instead, the cubic and force magnitude are
 * computed directly in each loop below:
 *   - Central force:  factor = GM*m / (r2 * r)   i.e. G*M*m/r^3
 *   - Pairwise force:  G*mi*mj / (dist2 * dist)  i.e. G*mi*mj/r^3
 * This avoids both the pow() libm call (~50-100 cycles) and the
 * function-call overhead on each of the ~8.4 billion invocations.
 */

void evolve(int count, double dt) {
  int step;
  int i, j;
  double size;
  double dt_over_mass;
  
  /* Precompute constant product to avoid repeated multiplication */
  const double GM_central = G * M_central;

  for(step = 1; step <= count; step++) {
    printf("timestep %d\n", step);
    printf("collisions %d\n", collisions);

    /* Compute viscosity and wind forces in a single fused pass.
     *
     * Original code called vis_forces() then wind_forces() in separate
     * loops per dimension (6 loops total), requiring multiple passes
     * over the force and viscosity arrays.
     *
     * Fused version: one loop per dimension computes both terms together.
     * f[j][i] = -vis[i] * (velo[j][i] + wind[j])
     * This halves the number of memory reads for the vis[] array. */
    for(j = 0; j < Ndim; j++) {
      double wj = wind[j];
      for(i = 0; i < Nbody; i++) {
        f[j][i] = -vis[i] * (velo[j][i] + wj);
      }
    }

    /* Fused central force: compute distance and apply force in one pass.
     *
     * Original code used 4 separate loops:
     *   1) Zero r[k]
     *   2) Accumulate norms via add_norms()
     *   3) Take sqrt()
     *   4) Apply central force via forces()
     *
     * Fusing into one loop keeps pos[] and r[] data in cache
     * and eliminates 3 extra passes over the arrays. */
    for(i = 0; i < Nbody; i++) {
      double r2 = pos[Xcoord][i] * pos[Xcoord][i] 
                + pos[Ycoord][i] * pos[Ycoord][i] 
                + pos[Zcoord][i] * pos[Zcoord][i];
      
      r[i] = sqrt(r2);
      
      /* r^3 = r^2 * r avoids calling pow() or extra sqrt() */
      double r3 = r2 * r[i];
      double factor = GM_central * mass[i] / r3;
      
      f[Xcoord][i] -= factor * pos[Xcoord][i];
      f[Ycoord][i] -= factor * pos[Ycoord][i];
      f[Zcoord][i] -= factor * pos[Zcoord][i];
    }

    /* Fused pairwise separation + force calculation.
     *
     * Original code stored ALL pairwise separations into delta_pos[]
     * and delta_r[] arrays first (Npair = ~8.4M entries), then looped
     * again to compute forces. This requires O(N^2) temporary storage
     * and causes severe cache thrashing for large N.
     *
     * Fused version computes separation, distance, and forces in one
     * pass. Local variables dx/dy/dz stay in registers, eliminating
     * all temporary array accesses. */
    for(i = 0; i < Nbody; i++) {
      /* Hoist invariant loads for particle i outside inner loop */
      double pxi = pos[Xcoord][i];
      double pyi = pos[Ycoord][i];
      double pzi = pos[Zcoord][i];
      double mi  = mass[i];
      double ri  = radius[i];
      double fxi = f[Xcoord][i];
      double fyi = f[Ycoord][i];
      double fzi = f[Zcoord][i];
      
      for(j = i + 1; j < Nbody; j++) {
        /* Compute separation vector */
        double dx = pxi - pos[Xcoord][j];
        double dy = pyi - pos[Ycoord][j];
        double dz = pzi - pos[Zcoord][j];
        
        /* Compute squared distance and distance */
        double dist2 = dx * dx + dy * dy + dz * dz;
        double dist = sqrt(dist2);
        
        size = ri + radius[j];
        
        /* Precompute force magnitude: G*m_i*m_j / r^3
         * Using dist2*dist = r^2 * r = r^3, avoids separate pow() call */
        double force_mag = G * mi * mass[j] / (dist2 * dist);
        
        /* Branchless direction: flip sign on collision.
         * Original used if/else which causes branch misprediction
         * (~15-20 cycle penalty) at ~50% rate in the inner loop.
         * Normal (attractive): direction = -1, so fxi += (-1)*F*dx => f[i] -= F
         * Collision (repulsive): direction = +1, so fxi += (+1)*F*dx => f[i] += F */
        double direction = -1.0 + 2.0 * (double)(dist < size);
        collisions += (dist < size);
        
        /* Apply pairwise forces using Newton's third law */
        double fx = direction * force_mag * dx;
        double fy = direction * force_mag * dy;
        double fz = direction * force_mag * dz;
        
        fxi += fx;
        fyi += fy;
        fzi += fz;
        f[Xcoord][j] -= fx;
        f[Ycoord][j] -= fy;
        f[Zcoord][j] -= fz;
      }
      /* Write back accumulated forces for particle i */
      f[Xcoord][i] = fxi;
      f[Ycoord][i] = fyi;
      f[Zcoord][i] = fzi;
    }

    /* Update positions using Euler integration */
    for(i = 0; i < Nbody; i++) {
      pos[Xcoord][i] += dt * velo[Xcoord][i];
      pos[Ycoord][i] += dt * velo[Ycoord][i];
      pos[Zcoord][i] += dt * velo[Zcoord][i];
    }

    /* Update velocities: precompute dt/mass[i] so that 3 divisions
     * per particle are replaced by 1 division + 3 multiplications.
     * Division costs ~20 cycles on Zen 2; multiplication ~4 cycles. */
    for(i = 0; i < Nbody; i++) {
      dt_over_mass = dt / mass[i];
      velo[Xcoord][i] += f[Xcoord][i] * dt_over_mass;
      velo[Ycoord][i] += f[Ycoord][i] * dt_over_mass;
      velo[Zcoord][i] += f[Zcoord][i] * dt_over_mass;
    }
  }
}
