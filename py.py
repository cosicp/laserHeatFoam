import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ── Parameters ────────────────────────────────────────────────────────────────
R        = 0.01    # cylinder radius [m]
z0       = 0.0     # starting height [m]
pitch    = 0  # helical pitch per revolution [m]
N        = 360     # number of laser points
n_arrows = 24      # how many normal arrows to draw

# ── Generate path ─────────────────────────────────────────────────────────────
revolutions = 2
thetas = np.linspace(0, revolutions * 2 * np.pi, N, endpoint=False)

x  =  R * np.cos(thetas)
y  =  R * np.sin(thetas)
z  =  z0 + pitch * (thetas / (2 * np.pi))

# Inward radial normals
nx = -np.cos(thetas)
ny = -np.sin(thetas)
nz =  np.zeros_like(thetas)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 8), facecolor='#0d0d0d')
ax  = fig.add_subplot(111, projection='3d', facecolor='#0d0d0d')

# Cylinder surface (wireframe guide)
theta_cyl = np.linspace(0, 2 * np.pi, 60)
z_cyl     = np.linspace(z.min() - 0.0002, z.max() + 0.0002, 2)
Tc, Zc    = np.meshgrid(theta_cyl, z_cyl)
Xc = R * np.cos(Tc)
Yc = R * np.sin(Tc)
ax.plot_surface(Xc, Yc, Zc, alpha=0.07, color='#2244aa', linewidth=0)
ax.plot_wireframe(Xc, Yc, Zc, alpha=0.12, color='#3355cc', linewidth=0.4, rstride=1, cstride=4)

# Laser path
ax.plot(x, y, z, color='#ffaa00', linewidth=1.8, label='Laser path', zorder=3)

# Normal vectors (subsample)
idx    = np.linspace(0, N - 1, n_arrows, dtype=int)
scale  = R * 0.45   # arrow length relative to radius

ax.quiver(
    x[idx], y[idx], z[idx],
    nx[idx] * scale, ny[idx] * scale, nz[idx] * scale,
    color='#00ffcc', linewidth=1.0, arrow_length_ratio=0.25,
    label='Inward normal  n=(−cosθ, −sinθ, 0)'
)

# Start / end markers
ax.scatter(*[x[0]],  *[y[0]],  *[z[0]],  color='#44ff44', s=60, zorder=5, label='Start')
ax.scatter(*[x[-1]], *[y[-1]], *[z[-1]], color='#ff4422', s=60, zorder=5, label='End')

# Axis labels & style
for spine in ['bottom', 'top', 'left', 'right']:
    pass  # 3D axes don't have spines the same way

ax.set_xlabel('x [m]', color='#aab', labelpad=8)
ax.set_ylabel('y [m]', color='#aab', labelpad=8)
ax.set_zlabel('z [m]', color='#aab', labelpad=8)
ax.tick_params(colors='#778', labelsize=8)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('#223')
ax.yaxis.pane.set_edgecolor('#223')
ax.zaxis.pane.set_edgecolor('#223')
ax.grid(True, color='#223344', linewidth=0.5)

ax.set_title('LPBF Radial Laser Path — Cylinder Interior',
             color='#7ec8ff', fontsize=13, pad=14)

legend = ax.legend(facecolor='#1a1a2e', edgecolor='#334', labelcolor='#ccd',
                   fontsize=9, loc='upper left')

# Inset: top-down view (XY plane) showing normal rotation
ax2 = fig.add_axes([0.72, 0.62, 0.22, 0.28], facecolor='#111122')
circle = plt.Circle((0, 0), R, color='#3355cc', fill=False, linewidth=1.2, alpha=0.5)
ax2.add_patch(circle)
ax2.plot(x, y, color='#ffaa00', linewidth=1.0, alpha=0.6)
ax2.quiver(
    x[idx], y[idx],
    nx[idx] * scale, ny[idx] * scale,
    color='#00ffcc', scale=1, scale_units='xy', angles='xy',
    width=0.008, headwidth=3, alpha=0.85
)
ax2.set_xlim(-R * 1.6, R * 1.6)
ax2.set_ylim(-R * 1.6, R * 1.6)
ax2.set_aspect('equal')
ax2.set_title('Top view (XY)', color='#7ec8ff', fontsize=8, pad=4)
ax2.tick_params(colors='#556', labelsize=6)
for sp in ax2.spines.values():
    sp.set_edgecolor('#334')

plt.tight_layout()
plt.savefig('cylinder_laser_path.png', dpi=150, bbox_inches='tight',
            facecolor=fig.get_facecolor())
plt.show()
print("Saved: cylinder_laser_path.png")