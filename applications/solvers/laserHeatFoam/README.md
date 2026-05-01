# laserHeatFoam

Transient heat conduction solver for Laser Powder Bed Fusion (LPBF) with
a volumetric laser source and progressive layer activation via
`fvMeshSubset`.

The base mesh contains all *N* pre-meshed layers; the solver activates
them one at a time, solving the heat equation only on the currently-active
sub-domain.

```
rho * cpEff * dT/dt - div(k * grad T) = Q
```

with `k = rp*kp + rm*km + rs*ks` reconstructed from the powder / melt /
re-solidified phase fractions, and `Q` from one of the supported laser
models (Gaussian surface, cylindrical, exponential, ellipsoidal, Goldak,
conical, Beer-Lambert, hybrid).

---

## 1. User guide

### Required case files

In addition to the standard `system/{controlDict,fvSchemes,fvSolution}`,
`constant/thermophysicalProperties`, and the laser inputs
(`constant/LaserProperties`, `constant/timeVsLaserPosition`,
`constant/timeVsLaserPower`), a multi-layer case needs:

#### `constant/multiLayerProperties`

```
baseplatePatch     baseplate;     // patch whose adjacent cells form layer 0
nLayers            4;             // number of layers to activate
cellsPerLayer      3;             // mesh cells through the thickness per layer
                                  // (optional, default 1)
layerDuration      200e-6;        // simulated time per layer [s]
layerInitialT      300;           // T applied to cells when activated [K]
exposedFacesPatch  topWall;       // optional: existing base-mesh patch onto
                                  // which the cut faces are placed. If
                                  // omitted, a synthetic 'oldInternalFaces'
                                  // patch is created (zeroGradient default).
```

#### Initial fields (`0/` or `initial/`)

`MUST_READ`: `T`, `rho`, `cp`, `Q`, `rc`. Reasonable starting values:

| Field | Typical |
|---|---|
| `T`   | uniform 300 K, `fixedValue 300` on the baseplate |
| `rho` | uniform value in the right ballpark; gets recomputed each step |
| `cp`  | uniform value; gets recomputed each step |
| `Q`   | uniform 0; rebuilt each step |
| `rc`  | 0 (powder) everywhere; rises monotonically as melting occurs |

### Mesh requirements

- The mesh must be **layered** in the build direction: face-neighbours
  of any layer's cells should be either in the same layer or in the
  layer immediately above. A mesh extruded normal to the baseplate
  patch (e.g. via `extrudeMesh` or a structured `blockMeshDict`)
  satisfies this.
- The total number of cells reachable through the thickness must be at
  least `nLayers * cellsPerLayer`. If not, the solver fails with a
  message naming the layer and step that ran out of cells.
- The mesh does **not** have to be axis-aligned in z; layer detection
  is purely topological.

### Time

`controlDict.endTime` is **overridden** by the solver. The actual end
time is `nLayers * layerDuration`. `deltaT`, `writeInterval`,
`adjustTimeStep`, and `maxDi` are honoured normally.

### Output

The solver writes the full base mesh. Cells in inactive layers retain
whatever was in the input field (typically `layerInitialT` or the
default values from `0/`); cells in the active sub-mesh carry the
running temperature solution. ParaView shows the sub-mesh growing layer
by layer if you threshold on the `meltHistory` or `rc` field.

The fields written each output step are: `T`, `rho`, `cp`, `rc`, `k`,
`cpEff`, `fL`, `rp`, `rm`, `rs`, `Q`, `meltHistory`, plus `gradT` (built
on demand at write-time).

### Worked example

`tutorials/multiLayer_basic` is a self-contained 200 × 800 × 120 µm
case with 12 layers, 3 cells per layer, a single Gaussian scan per
layer. Run with:

```bash
cd tutorials/multiLayer_basic
./Allclean
./Allrun
```

### Common issues

- **`Entry 'laplacian(k,subsetT)' not found`** — you wrote the solver
  yourself or modified the renaming logic; the sub-mesh field should be
  renamed back to its bare name (`T.rename("T")`) so existing
  `fvSchemes` entries apply.
- **`Layer k has no cells`** — the mesh has fewer reachable cell rows
  than `nLayers * cellsPerLayer`. Reduce `nLayers` or refine the mesh
  in the build direction.
- **Spurious heat sink at the top of the active region** — the synthetic
  `oldInternalFaces` patch defaults to a calculated/zero-gradient field;
  if that does not match your physics, define a real top patch on the
  base mesh with the BC you want and set `exposedFacesPatch` to it.

---

## 2. Developer guide

### Source layout

```
laserHeatFoam.C       Main: outer per-layer loop, sub-mesh shadowing,
                      base/sub field bookkeeping.
multiLayer.H          Reads multiLayerProperties, builds the per-layer
                      cell sets by face-neighbour growth.
createFields.H        Constructs all base-mesh fields (unchanged from
                      the single-layer ancestor).
updateThermo.H        T-dependent rho, cp, ks, km from polynomial fits.
updateLaser.H         Volumetric source Q from the configured laser
                      model and the current laser position/power.
updateCpeff.H         Apparent heat capacity (latent heat smearing).
updateFractions.H     rc, rp, rm, rs, meltHistory updates.
DiffusionNo.H         Diffusion-number CFL limiter (adjustTimeStep).
write.H               Per-step writer (gradT + runTime.write()).
```

### Layer detection (`multiLayer.H`)

Layers are constructed by face-neighbour expansion using
`fvMesh::cellCells()`:

1. **Seed for layer 0** = cells with a face on `baseplatePatch`,
   obtained from `polyPatch::faceCells()`.
2. **Seed for layer k > 0** = one halo step from `layerCells[k-1]`,
   excluding cells already in `covered`.
3. **Thicken to M cells** = `cellsPerLayer - 1` further halo steps,
   excluding `covered` and the partial layer being built.

The `haloStep` lambda factors out a single
"frontier → next, minus exclude sets" pass, used both for the seed
advance and for the thickening loop.

Empty halo steps abort the run with a message naming the layer and
step that failed, so configuration errors surface immediately.

### Per-layer scope (`laserHeatFoam.C`)

The trick that lets the existing `update*.H` includes work unchanged:
inside each layer iteration we open a nested scope that **C++-shadows**
`mesh` and the bare field names with sub-mesh equivalents.

```cpp
{
    const fvMesh& mesh = subMesher.subMesh();    // shadows base mesh

    volScalarField T(subMesher.interpolate(baseT));
    T.rename("T");                                // see 'naming' below
    // ... rho, cp, rc, meltHistory similarly ...

    // Derived fields (k, cpEff, Q, fL, fLold, rp, rm, rs, ks, km)
    // are constructed fresh on the sub-mesh.

    while (runTime.loop())
    {
        #include "updateThermo.H"
        #include "updateLaser.H"
        #include "DiffusionNo.H"
        #include "updateCpeff.H"
        #include "updateFractions.H"

        k = rp*kp + rm*km + rs*ks;

        while (simple.correctNonOrthogonal()) { TEqn.solve(); }

        // ... scatter sub-mesh values back to baseT, baseRho, ... ...

        {
            const fvMesh& mesh = baseMesh;        // un-shadow for write
            const volScalarField& T = baseT;
            #include "write.H"
        }
    }
}
```

Outside the scope, `T`, `mesh`, etc. resolve to the original base-mesh
objects re-bound as `baseT`, `baseMesh`, ...

### Field bookkeeping

| Category | Fields | Treatment |
|---|---|---|
| **Persistent** (carry state between layers) | `T`, `rho`, `cp`, `rc`, `meltHistory` | `subMesher.interpolate(...)` on entry, scatter via `cellMap()` after each `TEqn.solve()`. |
| **Derived** (rebuilt every timestep) | `k`, `cpEff`, `Q`, `fL`, `fLold`, `rp`, `rm`, `rs`, `ks`, `km` | Constructed fresh on the sub-mesh per layer. |

After every solve we scatter **every visualisable field** back to the
base mesh through `subMesher.cellMap()` so that `runTime.write()`
produces a complete picture of the build at each output time.

### Naming: the `subset*` rename

`fvMeshSubset::interpolate(vf)` returns a field with name
`"subset" + vf.name()`. Without renaming, `fvSchemes` entries keyed on
`laplacian(k,T)` would be looked up as `laplacian(k,subsetT)` and fail.
We call `subT.rename("T")` (etc.) right after construction so user
schemes / solvers continue to apply.

### Time handling

`runTime.setEndTime((layerI + 1) * layerDuration)` is called at the top
of each layer iteration. `runTime.loop()` then runs until that bumped
end time. Across layer boundaries, time keeps advancing — there is no
reset and no per-layer time directories.

Because the sub-mesh `T` is rebuilt fresh from the base field at each
layer transition, `T.oldTime()` for the previously-active cells equals
`T` on the first sub-step of a new layer (one step of `ddt` history is
dropped). This is acceptable for typical AM thermal simulations but is
worth bearing in mind if you investigate transients near layer
boundaries.

### Boundary on exposed faces

`fvMeshSubset::reset(activeCells, exposedPatchID)` controls where the
faces cut by the subset go:

- `exposedPatchID < 0` (default) → synthetic `oldInternalFaces` patch
  with calculated/zero-gradient-equivalent BCs.
- `exposedPatchID >= 0` → existing base-mesh patch; the BC carried by
  that patch (e.g. convective + radiative loss to ambient) is applied.

The choice is announced via `Info` so it is visible in the log.

### Things to watch when extending

- **Do not register sub-mesh fields with the same name on the same
  registry as base-mesh fields.** Field names live in the sub-mesh's
  `objectRegistry`; the base mesh has its own. Constructing a
  `volScalarField` named `T` on a sub-mesh does not collide with the
  base `T` because they use different registries.
- **Always scatter back before `runTime.write()`.** Otherwise the
  base-mesh fields written to disk will be stale by one timestep.
- **`simpleControl simple(baseMesh)` is fine outside the loop.**
  `correctNonOrthogonal()` only consults `system/fvSolution` via the
  case-level `objectRegistry`; the choice of mesh in the constructor
  does not affect the per-step iteration counter.

### Testing

Build:

```bash
wmake
```

Run the smoke test:

```bash
cd tutorials/multiLayer_basic
./Allclean && ./Allrun
grep -E "Activating layer|sub-mesh:" log.laserHeatFoam
```

Expected: one `Activating layer k` line per layer, sub-mesh cell count
growing monotonically by `(in-plane cells) * cellsPerLayer` each
activation, log ending in `End`.
