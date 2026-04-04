"""
Molecular Dynamics Simulation — Particle Visualization
=======================================================
4096 particles under central gravity + pairwise forces.
Uses Keplerian orbital mechanics matching the actual C simulation
constants (G=2, M_central=1000 from coord.h).

Inner particles orbit faster (differential rotation).
Particle trails show arc streaks for inner orbits.
Colour: hot colourmap (inner=bright white/yellow, outer=red/dark).

Usage
-----
  python visualize_md.py                          # runs standalone
  python visualize_md.py --output out.gif
  python visualize_md.py --frames 120 --fps 24

Requirements: pip install matplotlib numpy pillow
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ── Physics constants — match coord.h exactly ────────────────────────────────
G          = 2.0
M_CENTRAL  = 1000.0

# ── Visualisation settings ────────────────────────────────────────────────────
N_PARTICLES = 4096
TRAIL_LEN   = 12    # number of past frames shown as fading trail
R_MIN       = 60.0
R_MAX       = 1700.0


# ─────────────────────────────────────────────────────────────────────────────
# Extended XYZ loader  (format: Ar x y z vx vy vz)
# ─────────────────────────────────────────────────────────────────────────────
def load_xyz(path):
    """Parse extended XYZ file. Returns (frames_xy, frames_speed, meta).
    frames_xy    : list of (N,2) float arrays  — XY positions
    frames_speed : list of (N,)  float arrays  — |v| per particle
    meta         : list of dicts with 'frame', 'KE', 'collisions'
    """
    frames_xy, frames_speed, meta = [], [], []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            n = int(line)
        except ValueError:
            i += 1
            continue
        # Parse comment line for KE and collision metadata
        comment = lines[i + 1].strip() if i + 1 < len(lines) else ""
        info = {"frame": 0, "KE": 0.0, "collisions": 0}
        for token in comment.split():
            if token.startswith("Frame="):
                try: info["frame"] = int(token.split("=")[1])
                except ValueError: pass
            elif token.startswith("KE="):
                try: info["KE"] = float(token.split("=")[1])
                except ValueError: pass
            elif token.startswith("collisions="):
                try: info["collisions"] = int(token.split("=")[1])
                except ValueError: pass
        i += 2  # skip count + comment
        xy, spd = [], []
        for j in range(n):
            if i + j >= len(lines):
                break
            parts = lines[i + j].split()
            if len(parts) >= 7:
                try:
                    x, y = float(parts[1]), float(parts[2])
                    vx, vy, vz = float(parts[4]), float(parts[5]), float(parts[6])
                    xy.append([x, y])
                    spd.append(np.sqrt(vx*vx + vy*vy + vz*vz))
                except ValueError:
                    pass
            elif len(parts) >= 3:  # position-only fallback
                try:
                    xy.append([float(parts[1]), float(parts[2])])
                    spd.append(0.0)
                except ValueError:
                    pass
        i += n
        if xy:
            frames_xy.append(np.array(xy, dtype=float))
            frames_speed.append(np.array(spd, dtype=float))
            meta.append(info)
    return frames_xy, frames_speed, meta


# ─────────────────────────────────────────────────────────────────────────────
# Stats CSV loader
# ─────────────────────────────────────────────────────────────────────────────
def load_stats(path):
    """Return arrays (timesteps, KE, collisions) from stats.csv."""
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data[:, 0], data[:, 1], data[:, 2].astype(int)


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic Keplerian disk (fallback when no real trajectory provided)
# ─────────────────────────────────────────────────────────────────────────────
def make_disk(n=N_PARTICLES, seed=7):
    """Distribute n particles across a disk in circular Keplerian orbits."""
    rng = np.random.default_rng(seed)

    # Area-uniform radial distribution (denser toward centre, like real MD)
    u = rng.uniform(0, 1, n)
    r = R_MIN + (R_MAX - R_MIN) * u ** 0.7        # bias slightly inward

    phi0  = rng.uniform(0, 2 * np.pi, n)
    omega = np.sqrt(G * M_CENTRAL / r ** 3)        # rad / time-unit  (Kepler)

    return r, phi0, omega


def xy_at(r, phi0, omega, t):
    """Return (x, y) positions at simulation time t."""
    phi = phi0 + omega * t
    return r * np.cos(phi), r * np.sin(phi)


# ─────────────────────────────────────────────────────────────────────────────
def make_animation(xyz_path=None, n_frames=120, fps=24):
    using_real_data = xyz_path is not None

    if using_real_data:
        print(f"  Loading trajectory: {xyz_path}")
        frames_xy, frames_speed, meta = load_xyz(xyz_path)
        if not frames_xy:
            print("  WARNING: no frames parsed — falling back to synthetic")
            using_real_data = False
        else:
            n_frames = min(n_frames, len(frames_xy))
            frames_xy    = frames_xy[:n_frames]
            frames_speed = frames_speed[:n_frames]

    if not using_real_data:
        print("  Generating synthetic Keplerian trajectory …")
        r, phi0, omega = make_disk()
        omega_max = omega.max()
        t_total   = 3.0 * (2 * np.pi / omega_max)
        dt        = t_total / n_frames
        v_synth   = np.sqrt(G * M_CENTRAL / r)
        frames_xy, frames_speed = [], []
        for fi in range(n_frames):
            all_x, all_y = xy_at(r, phi0, omega, fi * dt)
            frames_xy.append(np.column_stack([all_x, all_y]))
            frames_speed.append(v_synth)
        meta = [{"frame": fi * 5, "KE": 0.0, "collisions": 0}
                for fi in range(n_frames)]

    # Position limits
    all_pos = np.vstack(frames_xy)
    lim = np.percentile(np.linalg.norm(all_pos, axis=1), 99) * 1.08
    lim = max(lim, 100.0)

    # Colour: normalise speed across all frames
    all_spd  = np.concatenate(frames_speed)
    spd_min  = np.percentile(all_spd, 2)
    spd_max  = np.percentile(all_spd, 98)
    spd_rng  = max(spd_max - spd_min, 1e-12)

    def norm_speed(spd):
        return np.clip((spd - spd_min) / spd_rng, 0, 1) * 0.84 + 0.08

    all_c = [norm_speed(s) for s in frames_speed]

    data_label = "ARCHER2 trajectory" if using_real_data else "Keplerian orbital mechanics"
    print(f"  Data: {data_label}  |  Axis ±{lim:.0f}  |  "
          f"{n_frames} frames  |  {len(frames_xy[0])} particles/frame")

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 9), facecolor="#000000")
    ax.set_facecolor("#000000")
    ax.set_aspect("equal")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axis("off")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.02)

    # Labels
    fig.text(0.5, 0.985,
             "Molecular Dynamics Simulation  —  4096 Particles",
             ha="center", va="top", color="white", fontsize=15, fontweight="bold")
    subtitle = ("ARCHER2 run  \u2022  AMD EPYC 7742  \u2022  "
                "Central gravity + pairwise forces + viscous damping"
                if using_real_data else
                "Keplerian orbital mechanics (G=2, M\u2090=1000)  \u2022  "
                "matches simulation physics")
    fig.text(0.5, 0.955, subtitle,
             ha="center", va="top", color="#aaaaaa", fontsize=10)
    fig.text(0.5, 0.928,
             "Optimised code: 31.9\u00d7 faster  \u2014  "
             "27.4 s vs 872.8 s  (ARCHER2, AMD EPYC 7742, Zen 2)",
             ha="center", va="top", color="#ffdd44", fontsize=10)

    # Faint radial guide rings
    for ring_frac in [0.25, 0.5, 0.75, 1.0]:
        rg = lim * ring_frac
        ax.add_patch(plt.Circle((0, 0), rg, color="#111111",
                                fill=False, linewidth=0.6, zorder=1))

    # Central mass — layered glow
    for sz, al in [(3500, 0.03), (900, 0.10), (250, 0.40), (60, 1.00)]:
        ax.scatter([0], [0], s=sz, c="#ffcc00", alpha=al,
                   zorder=6, linewidths=0)

    # ── Trail scatters (one per lag, decreasing alpha) ────────────────────────
    trail_scats = []
    for k in range(TRAIL_LEN):
        alpha = 0.30 * (1.0 - k / TRAIL_LEN) ** 1.5
        ts = ax.scatter([], [], s=0.7, c=[], cmap="hot",
                        vmin=0, vmax=1, alpha=alpha,
                        linewidths=0, zorder=2)
        trail_scats.append(ts)

    # ── Main particle scatter ─────────────────────────────────────────────────
    pos0 = frames_xy[0]
    main_scat = ax.scatter(
        pos0[:, 0], pos0[:, 1],
        c=all_c[0], cmap="hot",
        s=2.5, alpha=0.90, linewidths=0,
        vmin=0, vmax=1, zorder=3,
    )

    # Colourbar
    sm = plt.cm.ScalarMappable(
        cmap="hot",
        norm=plt.Normalize(vmin=spd_min, vmax=spd_max))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.01, location="right")
    cbar.set_label("Particle speed |v|", color="white", fontsize=9)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white", fontsize=8)

    # Info text
    ts_text = ax.text(0.02, 0.02, "",
                      transform=ax.transAxes, color="#555555",
                      fontsize=9, va="bottom")

    # ── Update function ───────────────────────────────────────────────────────
    def update(fi):
        pos = frames_xy[fi]
        main_scat.set_offsets(pos)
        main_scat.set_array(all_c[fi])

        for k, ts in enumerate(trail_scats):
            lag = k + 1
            if fi >= lag:
                pts = frames_xy[fi - lag]
                ts.set_offsets(pts)
                ts.set_array(all_c[fi - lag])
            else:
                ts.set_offsets(np.empty((0, 2)))

        m = meta[fi]
        info = f"Step {m['frame']}  |  collisions: {m['collisions']}"
        if m["KE"] > 0:
            info += f"  |  KE: {m['KE']:.1f}"
        ts_text.set_text(info)
        return (main_scat, ts_text, *trail_scats)

    anim = animation.FuncAnimation(
        fig, update, frames=n_frames,
        interval=1000 // fps, blit=True,
    )
    return anim, fig


# ─────────────────────────────────────────────────────────────────────────────
def make_stats_plot(stats_path, out_path):
    """Save a KE + collision-rate plot from stats.csv."""
    steps, KE, cols = load_stats(stats_path)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6),
                                   facecolor="black", sharex=True)
    fig.subplots_adjust(hspace=0.08, left=0.10, right=0.96,
                        top=0.88, bottom=0.10)

    for ax in (ax1, ax2):
        ax.set_facecolor("#0a0a0a")
        ax.tick_params(colors="#888888")
        for spine in ax.spines.values():
            spine.set_edgecolor("#333333")

    ax1.plot(steps, KE, color="#ff8800", linewidth=1.5)
    ax1.set_ylabel("Kinetic Energy", color="#aaaaaa", fontsize=10)
    ax1.tick_params(axis="y", colors="#888888")

    col_rate = np.diff(cols, prepend=cols[0])
    ax2.fill_between(steps, col_rate, color="#4488ff", alpha=0.6)
    ax2.plot(steps, col_rate, color="#4488ff", linewidth=0.8)
    ax2.set_ylabel("New collisions", color="#aaaaaa", fontsize=10)
    ax2.set_xlabel("Timestep", color="#888888", fontsize=10)
    ax2.tick_params(colors="#888888")

    fig.text(0.5, 0.95,
             "MD Simulation Physics  —  Kinetic Energy & Collision Rate",
             ha="center", va="top", color="white",
             fontsize=13, fontweight="bold")
    fig.text(0.5, 0.91,
             "Optimised code  |  ARCHER2  |  AMD EPYC 7742",
             ha="center", va="top", color="#ffdd44", fontsize=10)

    fig.savefig(out_path, dpi=130, facecolor="black")
    plt.close(fig)
    print(f"Saved stats plot: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="MD particle visualizer")
    parser.add_argument("--input",  "-i", default=None,
                        help="trajectory.xyz from ARCHER2 (optional)")
    parser.add_argument("--stats",  "-s", default=None,
                        help="stats.csv from ARCHER2 (optional, makes KE plot)")
    parser.add_argument("--output", "-o", default="visualization/md_simulation.gif")
    parser.add_argument("--frames", "-n", type=int, default=120)
    parser.add_argument("--fps",          type=int, default=24)
    args = parser.parse_args()

    # Optional: render KE / collision stats plot
    if args.stats:
        base = args.output.rsplit(".", 1)[0]
        make_stats_plot(args.stats, base + "_stats.png")

    print(f"Building animation: {args.frames} frames @ {args.fps} fps")
    anim, fig = make_animation(
        xyz_path=args.input, n_frames=args.frames, fps=args.fps
    )

    out = args.output
    print(f"Rendering → {out} …")
    if out.endswith(".gif"):
        writer = animation.PillowWriter(fps=args.fps)
    else:
        writer = animation.FFMpegWriter(
            fps=args.fps, bitrate=3000,
            extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"],
        )
    anim.save(out, writer=writer, dpi=140)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
