#!/bin/bash
#SBATCH --job-name=md_viz
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --account=

# -------------------------------------------------------
# Visualization data collection run
#
# Runs the optimised MD code for 10,000 timesteps to generate
# rich trajectory data for 3D portfolio visualization.
#
# Outputs (all in optimized/ directory after run):
#   trajectory.xyz  — extended XYZ: x y z vx vy vz per particle,
#                     200 frames sampled every 50 steps (~65 MB)
#   stats.csv       — timestep, KE, collision count per frame
#
# After job completes, download both files:
#   scp login.archer2.ac.uk:/path/to/optimized/trajectory.xyz visualization/
#   scp login.archer2.ac.uk:/path/to/optimized/stats.csv       visualization/
#
# Then render the visualization locally:
#   python visualization/visualize_md.py \
#       --input visualization/trajectory.xyz \
#       --stats visualization/stats.csv \
#       --output visualization/md_simulation.gif
# -------------------------------------------------------

module load PrgEnv-gnu

cd $SLURM_SUBMIT_DIR

echo "=========================================="
echo "MD Visualization Run"
echo "Start: $(date)"
echo "=========================================="

# Build
echo ""
echo "Building optimised code..."
cd optimized
make clean
make
echo "Build complete."

# Run — 10,000 timesteps (Nstep=1000, Nsave=10 hardcoded in control.c)
# Expected runtime: ~9 minutes on one ARCHER2 core
echo ""
echo "Running simulation: 10,000 timesteps (dt=0.02, total time=200)"
echo "Trajectory sampled every 50 steps -> 200 animation frames"
echo ""
time ./MD

# Report outputs
echo ""
echo "=========================================="
echo "Run complete: $(date)"
echo "=========================================="
echo ""
echo "Output files:"
ls -lh trajectory.xyz 2>/dev/null \
    && echo "  ✓ trajectory.xyz  (extended XYZ with velocities)" \
    || echo "  ✗ trajectory.xyz  NOT FOUND"
ls -lh stats.csv 2>/dev/null \
    && echo "  ✓ stats.csv       (KE + collision count per frame)" \
    || echo "  ✗ stats.csv       NOT FOUND"
echo ""
echo "Download commands (run on your laptop):"
echo "  scp \${USER}@login.archer2.ac.uk:\$(pwd)/trajectory.xyz visualization/"
echo "  scp \${USER}@login.archer2.ac.uk:\$(pwd)/stats.csv       visualization/"
