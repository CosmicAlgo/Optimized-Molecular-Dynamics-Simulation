# Molecular Dynamics Visualization

This directory contains tools for visualizing the MD simulation using Manim (Mathematical Animation Engine).

## Prerequisites

```bash
pip install -r requirements.txt
```

Or install Manim directly:
```bash
pip install manim
```

## Generating Trajectory Data

First, run the simulation with trajectory output enabled. The code has been modified to optionally save particle positions in XYZ format:

```bash
cd ../optimized
make clean && make
./MD
```

This will generate `trajectory.xyz` with particle positions sampled every 10 timesteps.

## Creating the Visualization

### Option 1: 3D Particle Simulation
```bash
manim -pqh md_visualization.py MDParticleSimulation
```

Flags:
- `-p`: Preview after rendering
- `-qh`: High quality (1080p)
- `-ql`: Low quality (faster, for testing)
- `-qk`: 4K quality

The output will be in `media/videos/md_visualization/720p30/` (or similar).

### Option 2: Speedup Comparison Chart
```bash
manim -pqh md_visualization.py MDSpeedupComparison
```

## Sample Data

A small sample trajectory is included in `sample_data/` for testing the visualization without running the full simulation.

## GitHub Integration

To embed the video in your README:

1. Upload the MP4 to GitHub (drag & drop into issue/PR or use GitHub Releases)
2. Or convert to GIF using:
   ```bash
   ffmpeg -i MDParticleSimulation.mp4 -vf "fps=30,scale=480:-1:flags=lanczos,split[s0][s1];[s0]palettegen=[s1]paletteuse=dither=bayer" -loop 0 output.gif
   ```
3. Reference in README:
   ```markdown
   ![MD Simulation](media/videos/md_visualization/720p30/MDParticleSimulation.mp4)
   ```

## Performance Notes

- Full 4096-particle animation may be slow. The script limits to 500 particles for smooth playback.
- Render time: ~5-10 minutes for high quality on modern hardware.
- For quick tests, use `-ql` (low quality) flag.
