---
name: laserHeatFoam
description: Use this skill when working in this directory. This project is an OpenFOAM finite-volume solver named laserHeatFoam, with solver sources, header includes, and Make files following OpenFOAM application conventions.
---

# laserHeatFoam OpenFOAM Solver

This directory is an OpenFOAM solver application, not a standalone C++ project.

## Project Context

- Solver executable/source: `laserHeatFoam.C`
- OpenFOAM include fragments: `createFields.H`, `updateThermo.H`, `updateLaser.H`, `updateLaserCyl.H`, `updateLaserExp.H`, `updateConvectionVol.H`, `updateCpeff.H`, `write.H`, `DiffusionNo.H`
- Build configuration: `Make/files` and `Make/options`
- Expected build tool: `wmake`, from a shell where OpenFOAM has been sourced.

## Working Rules

- Preserve OpenFOAM coding patterns already present in the solver.
- Treat `.H` files here as source fragments included by `laserHeatFoam.C`, not normal public C++ headers.
- Prefer OpenFOAM field, dimensioned type, mesh, dictionary, and finite-volume APIs over generic C++ replacements.
- Keep changes scoped to solver behavior unless the user asks for broader library or case changes.
- Do not assume ordinary CMake, Makefile, or package-manager workflows apply.

## Verification

When possible, verify edits with:

```sh
wmake
```

If OpenFOAM is not sourced, report that `wmake` cannot run in the current shell and leave the code changes ready for an OpenFOAM environment.
