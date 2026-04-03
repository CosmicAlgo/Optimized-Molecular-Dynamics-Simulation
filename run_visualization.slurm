#!/bin/bash
#SBATCH --job-name=md_viz
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --account=

# -------------------------------------------------------
# Visualization run: generates trajectory.xyz for Manim
# 
# This script runs the optimized MD code to generate
# particle trajectory data for 3D visualization.
# 
# Output: trajectory.xyz (XYZ format, sampled every 5 timesteps)
# 
# After completion:
#   1. Download trajectory.xyz to your laptop
#   2. Run: manim -pqh visualization/md_visualization.py MDParticleSimulation
# -------------------------------------------------------

module load PrgEnv-gnu

# Change to submission directory
cd $SLURM_SUBMIT_DIR

echo "=========================================="
echo "MD Visualization Run"
echo "Start: $(date)"
echo "=========================================="

# Build optimized code
echo "Building optimized code..."
cd optimized
make clean
make

# Run simulation (500 timesteps = full run)
# This generates trajectory.xyz
echo ""
echo "Running simulation (500 timesteps)..."
echo "Trajectory sampled every 5 timesteps"
echo ""

./MD

# Check output
echo ""
echo "=========================================="
echo "Run complete: $(date)"
echo "=========================================="
echo ""
echo "Output files:"
ls -lh trajectory.xyz 2>/dev/null && echo "  ✓ trajectory.xyz ready for download" || echo "  ✗ trajectory.xyz not found"
echo ""
echo "Next steps:"
echo "  1. Download trajectory.xyz to your laptop"
echo "  2. Install Manim: pip install manim"
echo "  3. Generate video: manim -pqh visualization/md_visualization.py MDParticleSimulation"
