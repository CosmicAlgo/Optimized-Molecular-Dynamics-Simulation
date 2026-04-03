"""
Molecular Dynamics 3D Visualization using Manim

This script creates a professional 3D animation of the MD simulation.
Requires: manim, numpy

Usage:
    manim -pqh md_visualization.py MDParticleSimulation
    
For lower quality (faster):
    manim -pql md_visualization.py MDParticleSimulation
"""

from manim import *
import numpy as np


class MDParticleSimulation(ThreeDScene):
    """
    3D visualization of molecular dynamics simulation.
    Shows 4096 particles evolving under gravitational forces.
    """
    
    def construct(self):
        # Set up 3D camera
        self.set_camera_orientation(phi=75 * DEGREES, theta=-45 * DEGREES)
        
        # Title and info
        title = Text(
            "Molecular Dynamics Simulation\n4096 particles | 500 timesteps",
            font_size=24,
            color=WHITE
        ).to_corner(UL)
        
        stats = Text(
            "31.9x speedup on ARCHER2",
            font_size=20,
            color=YELLOW
        ).to_corner(UR)
        
        self.add_fixed_in_frame_mobjects(title, stats)
        
        # Create coordinate axes
        axes = ThreeDAxes(
            x_range=[-2, 2, 0.5],
            y_range=[-2, 2, 0.5],
            z_range=[-2, 2, 0.5],
            x_length=4,
            y_length=4,
            z_length=4,
        )
        
        axes_label = Text("Simulation Space (arbitrary units)", font_size=16)
        axes_label.next_to(axes, DOWN)
        
        self.add_fixed_in_frame_mobjects(axes_label)
        
        # Load trajectory data
        # Format: timestep, x, y, z, vx, vy, vz (one line per particle per frame)
        try:
            trajectory = self.load_trajectory("trajectory.xyz")
        except FileNotFoundError:
            # Generate synthetic data for demo if trajectory not available
            trajectory = self.generate_synthetic_trajectory()
        
        num_frames = len(trajectory)
        num_particles = len(trajectory[0])
        
        # Create particle dots
        # Use small dots for performance with 4096 particles
        particles = VGroup(*[
            Dot3D(
                point=[0, 0, 0],
                radius=0.02,
                color=BLUE if i < num_particles // 2 else RED,
                glow=False
            )
            for i in range(min(num_particles, 500))  # Limit for animation smoothness
        ])
        
        self.add(axes, particles)
        
        # Add central mass representation
        central_mass = Sphere(
            center=[0, 0, 0],
            radius=0.05,
            color=YELLOW,
            gloss=0.8
        )
        self.add(central_mass)
        
        # Animation
        self.begin_ambient_camera_rotation(rate=0.1)
        
        # Animate through frames
        frame_duration = 0.1  # seconds per frame
        
        for frame_idx in range(0, min(num_frames, 100), 2):  # Skip frames for speed
            frame_data = trajectory[frame_idx]
            
            # Update particle positions
            new_positions = []
            for i, dot in enumerate(particles):
                if i < len(frame_data):
                    x, y, z = frame_data[i][:3]
                    # Scale for visualization
                    new_pos = [x * 2, y * 2, z * 2]
                    new_positions.append(new_pos)
                else:
                    new_positions.append(dot.get_center())
            
            # Animate to new positions
            self.play(
                *[dot.animate.move_to(pos) for dot, pos in zip(particles, new_positions)],
                run_time=frame_duration,
                rate_func=linear
            )
        
        self.wait(2)
        self.stop_ambient_camera_rotation()
    
    def load_trajectory(self, filename):
        """Load trajectory from XYZ file format."""
        frames = []
        current_frame = []
        
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            # Read number of atoms
            if lines[i].strip().isdigit():
                n_atoms = int(lines[i].strip())
                i += 1  # Skip comment line
                
                # Read atom positions
                frame = []
                for j in range(n_atoms):
                    if i + j < len(lines):
                        parts = lines[i + j].split()
                        if len(parts) >= 4:
                            # Format: element x y z
                            x, y, z = map(float, parts[1:4])
                            frame.append([x, y, z])
                
                frames.append(frame)
                i += n_atoms
            else:
                i += 1
        
        return frames
    
    def generate_synthetic_trajectory(self, num_frames=50, num_particles=500):
        """Generate synthetic orbital trajectory for demo purposes."""
        frames = []
        
        for t in range(num_frames):
            frame = []
            for i in range(num_particles):
                # Create orbital motion around center
                angle = 2 * np.pi * i / num_particles + t * 0.1
                radius = 0.5 + 0.3 * np.sin(2 * np.pi * i / num_particles)
                
                x = radius * np.cos(angle)
                y = radius * np.sin(angle)
                z = 0.1 * np.sin(4 * angle + t * 0.2)
                
                frame.append([x, y, z])
            
            frames.append(frame)
        
        return frames


class MDSpeedupComparison(Scene):
    """
    2D animation showing speedup comparison between original and optimized.
    """
    
    def construct(self):
        # Title
        title = Text("Performance Optimization Results", font_size=32)
        title.to_edge(UP)
        
        # Speedup data
        data = {
            "Original (O0)": 872.8,
            "Optimized (O0)": 108.3,
            "Original (O3)": 68.3,
            "Optimized (O3)": 27.4
        }
        
        # Create bar chart
        colors = [RED, ORANGE, BLUE, GREEN]
        bars = VGroup()
        labels = VGroup()
        
        max_time = max(data.values())
        scale = 5 / max_time
        
        for i, (name, time) in enumerate(data.items()):
            # Bar
            bar = Rectangle(
                height=time * scale,
                width=1.2,
                fill_color=colors[i],
                fill_opacity=0.8,
                stroke_color=WHITE,
                stroke_width=2
            )
            bar.shift(DOWN * 2.5 + RIGHT * (i * 2 - 3) + UP * time * scale / 2)
            
            # Label
            label = Text(name, font_size=16)
            label.next_to(bar, DOWN, buff=0.2)
            
            # Time value
            time_text = Text(f"{time:.1f}s", font_size=18)
            time_text.next_to(bar, UP, buff=0.2)
            
            bars.add(bar)
            labels.add(label)
            labels.add(time_text)
        
        # Speedup annotation
        speedup = Text("31.9x Speedup", font_size=36, color=YELLOW)
        speedup.to_edge(DOWN)
        
        # Animate
        self.play(Write(title))
        self.play(Create(bars), run_time=2)
        self.play(Write(labels))
        self.play(Write(speedup))
        
        self.wait(2)
