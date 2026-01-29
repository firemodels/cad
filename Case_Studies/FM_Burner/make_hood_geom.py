#!/usr/bin/env python3
"""
make_hood_geom.py

Closed, watertight, manifold hood solid with thickness everywhere.

Key topology (matches your note):
  - Inner cylinder (duct inner wall) STOPS at ZMAX (cannot overshoot cap bottom).
  - Outer cylinder CONTINUES to ZMAX+CAP_THK (cap top is part of outer cylinder).
  - Bottom gas-facing cap is ONLY the disk 0..Ri_top at ZMAX (normal -z).
      * This disk gets SURF_ID(2) = HOOD_OUTFLOW.
  - Cap top face is a full disk 0..Ro_top at ZMAX+CAP_THK (normal +z), INERT.
  - No coincident opposite-facing faces (avoids non-manifold edges).

Output:
  hood.txt
"""

from __future__ import annotations
import numpy as np

# -----------------------------
# Parameters
# -----------------------------
NTHETA = 48  # a bit finer helps reduce skinny triangles

# Radii (m)
Ri_bot = 0.3048          # 24 in ID / 2
Ri_top = 0.0762          #  6 in ID / 2
t      = 0.02            # 2 cm thickness
Ro_bot = Ri_bot + t
Ro_top = Ri_top + t

# Z locations (m)
Z_bot          = 0.8128
Z_frustum_top  = Z_bot + 0.3048   # 12 in transition (adjust if desired)
ZMAX           = 1.40             # <-- SET THIS to your chamber zmax
CAP_THK        = 0.02             # 2 cm cap thickness above ZMAX
Z_cap_top      = ZMAX + CAP_THK

OUTFILE = "hood.txt"

# SURF indices
SURF_INERT = 1
SURF_OUT   = 2

# -----------------------------
# Geometry helpers
# -----------------------------
def ring(r: float, z: float) -> np.ndarray:
    th = np.linspace(0.0, 2.0*np.pi, NTHETA, endpoint=False)
    return np.column_stack((r*np.cos(th), r*np.sin(th), np.full(NTHETA, z)))

verts: list[np.ndarray] = []
faces: list[tuple[int,int,int,int]] = []

def add_ring(r: float, z: float) -> np.ndarray:
    i0 = len(verts)
    verts.extend([v for v in ring(r, z)])
    return np.arange(i0, i0+NTHETA, dtype=int)

def add_vertex(x: float, y: float, z: float) -> int:
    i = len(verts)
    verts.append(np.array([x, y, z], float))
    return i

def wrap(i: int) -> int:
    return i % NTHETA

def add_tri(a: int, b: int, c: int, surf: int) -> None:
    faces.append((a, b, c, surf))

def add_quad(a: int, b: int, c: int, d: int, surf: int) -> None:
    # two triangles
    add_tri(a, b, c, surf)
    add_tri(a, c, d, surf)

def force_outward(a: int, b: int, c: int, desired: np.ndarray) -> tuple[int,int,int]:
    va, vb, vc = verts[a], verts[b], verts[c]
    n = np.cross(vb - va, vc - va)
    return (a, c, b) if np.dot(n, desired) < 0.0 else (a, b, c)

def radial_dir(a: int, b: int, c: int) -> np.ndarray:
    p = (verts[a] + verts[b] + verts[c]) / 3.0
    r = np.array([p[0], p[1], 0.0])
    nr = np.linalg.norm(r)
    if nr < 1e-12:
        return np.array([1.0, 0.0, 0.0])
    return r / nr

def add_tri_out(a: int, b: int, c: int, surf: int, desired: np.ndarray) -> None:
    aa, bb, cc = force_outward(a, b, c, desired)
    add_tri(aa, bb, cc, surf)

def add_quad_out(a: int, b: int, c: int, d: int, surf: int, desired: np.ndarray) -> None:
    add_tri_out(a, b, c, surf, desired)
    add_tri_out(a, c, d, surf, desired)

# -----------------------------
# Rings
# -----------------------------
ri_bot = add_ring(Ri_bot, Z_bot)
ro_bot = add_ring(Ro_bot, Z_bot)

ri_mid = add_ring(Ri_top, Z_frustum_top)
ro_mid = add_ring(Ro_top, Z_frustum_top)

# duct at ZMAX
ri_top = add_ring(Ri_top, ZMAX)
ro_top = add_ring(Ro_top, ZMAX)

# outer ring at cap top (outer cylinder continues)
ro_cap = add_ring(Ro_top, Z_cap_top)

# centers for disks
c_out = add_vertex(0.0, 0.0, ZMAX)       # outflow disk center (bottom of cap)
c_top = add_vertex(0.0, 0.0, Z_cap_top)  # top cap disk center

# -----------------------------
# Frustum walls
# -----------------------------
for i in range(NTHETA):
    ip = wrap(i+1)

    # Outer frustum (+r)
    desired = radial_dir(ro_bot[i], ro_mid[i], ro_mid[ip])
    add_quad_out(ro_bot[i], ro_mid[i], ro_mid[ip], ro_bot[ip], SURF_INERT, desired)

    # Inner frustum (-r)
    desired = -radial_dir(ri_bot[i], ri_mid[i], ri_mid[ip])
    add_quad_out(ri_bot[i], ri_bot[ip], ri_mid[ip], ri_mid[i], SURF_INERT, desired)

# -----------------------------
# Duct walls up to ZMAX
# -----------------------------
for i in range(NTHETA):
    ip = wrap(i+1)

    # Outer duct (+r) from ro_mid -> ro_top
    desired = radial_dir(ro_mid[i], ro_top[i], ro_top[ip])
    add_quad_out(ro_mid[i], ro_top[i], ro_top[ip], ro_mid[ip], SURF_INERT, desired)

    # Inner duct (-r) from ri_mid -> ri_top (STOP at ZMAX)
    desired = -radial_dir(ri_mid[i], ri_top[i], ri_top[ip])
    add_quad_out(ri_mid[i], ri_mid[ip], ri_top[ip], ri_top[i], SURF_INERT, desired)

# -----------------------------
# Bottom annulus closure at Z_bot (normal -z)
# -----------------------------
mz = np.array([0.0, 0.0, -1.0])
for i in range(NTHETA):
    ip = wrap(i+1)
    add_quad_out(ri_bot[i], ri_bot[ip], ro_bot[ip], ro_bot[i], SURF_INERT, mz)

# -----------------------------
# Outer cylinder continues upward: ro_top -> ro_cap (+r)
# (This is the "cap is part of outer cylinder" requirement.)
# -----------------------------
for i in range(NTHETA):
    ip = wrap(i+1)
    desired = radial_dir(ro_top[i], ro_cap[i], ro_cap[ip])
    add_quad_out(ro_top[i], ro_cap[i], ro_cap[ip], ro_top[ip], SURF_INERT, desired)

# -----------------------------
# Top cap face at Z_cap_top (normal +z), full disk 0..Ro_top
# -----------------------------
pz = np.array([0.0, 0.0, +1.0])
for i in range(NTHETA):
    ip = wrap(i+1)
    add_tri_out(c_top, ro_cap[i], ro_cap[ip], SURF_INERT, pz)

# -----------------------------
# Bottom of cap at ZMAX: ONLY the disk 0..Ri_top (normal -z)
# This connects to inner cylinder and avoids non-manifold edges.
# -----------------------------
for i in range(NTHETA):
    ip = wrap(i+1)
    # triangles should point downward (-z) as outward normal of solid (solid is above)
    add_tri_out(c_out, ri_top[ip], ri_top[i], SURF_OUT, mz)

# -----------------------------
# Write GEOM
# -----------------------------
with open(OUTFILE, "w") as f:
    f.write("&GEOM ID='HOOD',CELL_BLOCK_IOR=3,\n")
    f.write("      SURF_ID='INERT','HOOD_OUTFLOW',\n")
    f.write("      VERTS=\n")
    for v in verts:
        f.write("      %.6f, %.6f, %.6f,\n" % (v[0], v[1], v[2]))
    f.write("      FACES=\n")
    for a,b,c,s in faces:
        f.write("      %d, %d, %d, %d,\n" % (a+1, b+1, c+1, s))
    f.write(" /\n")

print(f"Wrote {OUTFILE}")

