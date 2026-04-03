# Molecular Dynamics Simulation — HPC Performance Optimisation

A case study in performance engineering applied to an $N$-body molecular dynamics simulation ($N = 4096$ particles, 500 timesteps) running on a single core of the [ARCHER2](https://www.archer2.ac.uk/) supercomputer (AMD EPYC 7742, Zen 2).

**Result: 31.9× speedup** over the unoptimised baseline (874.8 s → 27.4 s), achieved through seven source-level transformations combined with targeted compiler flags — no parallelism used.

---

## Structure

```
original/        Original unmodified simulation code
optimized/       Optimised version + SLURM benchmark scripts
Test/            diff-output correctness checker
Coursework2/     Benchmark results and LaTeX report
input.dat        Simulation input (4096 particles)
```

---

## Optimisations Applied

| # | Transformation | Key Benefit |
|---|---------------|-------------|
| 1 | Replace `pow(rv, 3.0)` with `r*r*r` + inline | Eliminates ~50–100 cycle libm call per pair |
| 2 | Fuse viscosity + wind force loops | Halves reads of `vis[]`, removes 3000 function calls/timestep |
| 3 | Fuse central force (4 loops → 1) | Keeps `pos[]`/`mass[]` in L1 cache |
| 4 | Fuse pairwise separation + force calc | Eliminates ~256 MB temporary arrays, removes cache thrashing |
| 5 | Branchless collision detection | Removes ~50% branch mispredictions in inner loop |
| 6 | Strength reduction in velocity update | 3 divisions → 1 division + 3 multiplications per particle |
| 7 | Remove unused temporary array allocation | Saves ~256 MB heap, reduces page faults and TLB pressure |

The largest single gain comes from **Opt 4** (fused pairwise loop), which eliminates the need to store all $N_\text{pair} = 8{,}386{,}560$ pairwise separation vectors before computing forces — a design that caused the working set to overflow the 16 MB L3 cache on every timestep.

---

## Benchmark Results (ARCHER2, single core)

| Source | Compiler flags | Time (s) | Speedup |
|--------|---------------|----------|---------|
| Original | `-O0 -g` | 874.8 ± 2.8 | 1.0× |
| Optimised | `-O0 -g` | 87.1 ± 0.3 | 10.0× |
| Original | production | 68.3 ± 0.8 | 12.8× |
| **Optimised** | **production** | **27.4 ± 0.1** | **31.9×** |

Production flags: `-O3 -march=znver2 -ffast-math -funroll-loops -flto`  
Compiler: GCC 11.2.0 (HPE, via `PrgEnv-gnu`)  
Each entry is the mean ± std dev of 5 independent runs.

The source-level changes and compiler flags interact non-additively ($10.0× \times 12.8× \neq 31.9×$): changes like strength reduction overlap with what `-O3` already does, while eliminating the temporary arrays is beyond the compiler's reach.

---

## Profiling (gprof, baseline `-O0 -g -pg`)

| Function | % Time | Calls |
|----------|--------|-------|
| `evolve` | 48.3% | 5 |
| `forces` | 41.3% | ~8 × 10⁹ |
| `add_norms` | 10.5% | 3,000 |
| `vis_forces` / `wind_forces` | < 0.1% | 1,500 each |

`forces()` body: `return Wv * deltav / pow(rv, 3.0)` — a single libm call dominating 41% of runtime across 8 billion invocations.

---

## Running the Benchmarks (ARCHER2)

```bash
# Set your account in the SLURM scripts first:
# sed -i 's/YOUR_ACCOUNT/your-archer2-account/' optimized/*.slurm

# Baseline (original code, -O0): ~15 min
sbatch optimized/run_baseline.slurm

# Optimised (production flags): ~5 min
sbatch optimized/run_optimized.slurm

# Full comparison (all variants + gprof + correctness): ~6 hours
sbatch optimized/run_comparison.slurm

# Quick correctness verification only: ~5 min
sbatch optimized/run_verify.slurm
```

Results are written to `results/` in the submission directory.

---

## Correctness

The `Test/diff-output` utility computes per-field relative error between simulation outputs and flags any value exceeding 5%. The optimised code matches the original to within floating-point reassociation tolerance introduced by `-ffast-math` (well under 0.1%).

---

## Report

The full write-up (methodology, profiling analysis, per-optimisation discussion, results, code quality) is in [`Coursework2/report_cw2/report.tex`](Coursework2/report_cw2/report.tex).

---

## Build

```bash
# Optimised build (production flags)
cd optimized && make

# Original build
cd original && make

# Correctness checker
cd Test && make
```

---

## 3D Visualization (Manim)

A professional 3D animation of the simulation can be generated using [Manim](https://www.manim.community/) (Mathematical Animation Engine).

### Quick Start

```bash
# Install Manim
pip install -r visualization/requirements.txt

# Run simulation to generate trajectory data
cd optimized
make && ./MD
# Creates trajectory.xyz with particle positions

# Generate 3D animation
manim -pqh visualization/md_visualization.py MDParticleSimulation
# Output: media/videos/md_visualization/720p30/MDParticleSimulation.mp4
```

### Features

- **3D particle simulation** with 4096 particles in orbital motion
- **Real-time camera rotation** showing spatial dynamics
- **Performance stats overlay** (31.9× speedup highlighted)
- **Speedup comparison chart** (alternative 2D scene)

See [`visualization/README.md`](visualization/README.md) for detailed usage and embedding in GitHub.

---

## Portfolio

This project demonstrates:
- **Performance analysis**: profiling with gprof, cache-aware optimization
- **Low-level optimization**: branchless programming, strength reduction, loop fusion
- **HPC systems**: ARCHER2, SLURM, compiler toolchains
- **Visualization**: 3D scientific animation with Manim
