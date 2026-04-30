/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laserHeatFoam

Group
    LPBF solvers

Description
    Heat conduction solver for Laser Powder Bed Fusion (LPBF) with a
    volumetric heat source and progressive layer activation via
    fvMeshSubset. The base mesh contains all N layers; at each outer
    iteration the active cell set grows by one layer, a sub-mesh is
    rebuilt, the transient heat equation is solved on the sub-mesh for
    layerDuration seconds, and the persistent state is scattered back
    to the base mesh for output.

    Layer activation is topological: layer 0 is the set of cells touching
    a user-named baseplate patch; layer k is the set of face-neighbours
    of layer k-1 that have not yet been activated. See multiLayer.H.

    \heading Solver details
    The solver solves the transient heat conduction equation:

    \f[
        \rho c_{p}^{eff} \ddt{T}  = \div \left( k \grad T \right) + Q
    \f]

Author
    Petar Cosic, UCD.
    Multi-layer extension by Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "interpolationTable.H"
#include "mathematicalConstants.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Heat conduction solver for Laser Powder Bed Fusion (LPBF) "
        "with a volumetric heat source and progressive layer activation."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "multiLayer.H"

    // Re-bind base-mesh fields under explicit names so we can shadow the
    // bare names with sub-mesh equivalents inside the layer loop without
    // losing access to the base versions.
    const fvMesh& baseMesh = mesh;
    volScalarField& baseT = T;
    volScalarField& baseRho = rho;
    volScalarField& baseCp = cp;
    volScalarField& baseRc = rc;
    volScalarField& baseK = k;
    volScalarField& baseCpEff = cpEff;
    volScalarField& baseFL = fL;
    volScalarField& baseRp = rp;
    volScalarField& baseRm = rm;
    volScalarField& baseRs = rs;
    volScalarField& baseQ = Q;
    volScalarField& baseMeltHistory = meltHistory;

    fvMeshSubset subMesher(baseMesh);
    labelHashSet activeCells(baseMesh.nCells());

    Info<< nl
        << "Simulating " << nLayers << " layer(s); each layer runs for "
        << layerDuration << " s." << nl
        << "Total simulated time = " << nLayers*layerDuration << " s." << nl
        << "(controlDict endTime is overridden per-layer.)" << nl << endl;

    Info<< "Calculating temperature distribution" << nl << endl;

    for (label layerI = 0; layerI < nLayers; ++layerI)
    {
        Info<< "Activating layer " << layerI
            << " (" << layerCells[layerI].size() << " new cells)" << endl;

        // Initialise temperature on cells being activated this iteration.
        // Already-active cells keep their running thermal history.
        for (const label cellI : layerCells[layerI])
        {
            baseT[cellI] = layerInitialT;
        }
        baseT.correctBoundaryConditions();

        // Grow the active set and rebuild the sub-mesh.
        activeCells |= layerCells[layerI];
        subMesher.reset(activeCells, exposedPatchID);

        // Extend simulation end time to cover this layer.
        runTime.setEndTime((layerI + 1)*layerDuration);

        Info<< "    sub-mesh: " << subMesher.subMesh().nCells()
            << " cells; layer end time: " << runTime.endTime().value()
            << " s" << nl << endl;

        // Sub-mesh scope: the bare names mesh, T, rho, ... are shadowed
        // here so that the existing update*.H, DiffusionNo.H and write.H
        // includes operate on the sub-mesh without modification.
        {
            const fvMesh& mesh = subMesher.subMesh();

            // Persistent state copied from the base mesh.
            // fvMeshSubset::interpolate names the result "subset<base>";
            // rename so fvSchemes/fvSolution entries keyed on the bare
            // field name (e.g. laplacian(k,T)) continue to apply.
            volScalarField T(subMesher.interpolate(baseT));
            T.rename("T");
            volScalarField rho(subMesher.interpolate(baseRho));
            rho.rename("rho");
            volScalarField cp(subMesher.interpolate(baseCp));
            cp.rename("cp");
            volScalarField rc(subMesher.interpolate(baseRc));
            rc.rename("rc");
            volScalarField meltHistory
            (
                subMesher.interpolate(baseMeltHistory)
            );
            meltHistory.rename("meltHistory");

            // Derived fields rebuilt every step from T, rc, ...; constructed
            // fresh on the sub-mesh here. Dimensions and IO flags mirror
            // createFields.H.
            volScalarField k
            (
                IOobject
                (
                    "k",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "k",
                    dimensionSet(1, 1, -3, -1, 0, 0, 0),
                    0.0
                )
            );

            volScalarField cpEff
            (
                IOobject
                (
                    "cpEff",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                cp
            );

            volScalarField Q
            (
                IOobject
                (
                    "Q",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "Q",
                    dimensionSet(1, -1, -3, 0, 0, 0, 0),
                    0.0
                )
            );

            volScalarField fL
            (
                IOobject
                (
                    "fL",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            );

            volScalarField fLold
            (
                IOobject
                (
                    "fLold",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fL
            );

            volScalarField rp
            (
                IOobject
                (
                    "rp",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            );

            volScalarField rm
            (
                IOobject
                (
                    "rm",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            );

            volScalarField rs
            (
                IOobject
                (
                    "rs",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            );

            volScalarField ks
            (
                IOobject
                (
                    "ks",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "ks",
                    dimensionSet(1, 1, -3, -1, 0, 0, 0),
                    0.0
                )
            );

            volScalarField km
            (
                IOobject
                (
                    "km",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "km",
                    dimensionSet(1, 1, -3, -1, 0, 0, 0),
                    0.0
                )
            );

            // Time loop for this layer.
            while (runTime.loop())
            {
                Info<< "Time = " << runTime.timeName() << nl << endl;

                // Preserve the previous liquid fraction before thermo
                // updates.
                fLold = fL;

                #include "updateThermo.H"
                #include "updateLaser.H"
                #include "DiffusionNo.H"
                #include "updateCpeff.H"
                #include "updateFractions.H"

                // Effective conductivity from the current phase mix.
                k = rp*kp + rm*km + rs*ks;

                // Solve the transient heat equation with an apparent heat
                // capacity.
                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix TEqn
                    (
                        rho*cpEff*fvm::ddt(T) - fvm::laplacian(k, T) == Q
                    );

                    TEqn.solve();
                }

                Info<< "T (active sub-mesh): max = " << max(T).value()
                    << ", min = " << min(T).value() << endl;

                // Scatter sub-mesh state back to the base mesh so that
                // write.H produces complete output for the whole domain.
                {
                    const labelList& cMap = subMesher.cellMap();
                    forAll(cMap, i)
                    {
                        const label bI = cMap[i];

                        baseT[bI] = T[i];
                        baseRho[bI] = rho[i];
                        baseCp[bI] = cp[i];
                        baseRc[bI] = rc[i];
                        baseK[bI] = k[i];
                        baseCpEff[bI] = cpEff[i];
                        baseFL[bI] = fL[i];
                        baseRp[bI] = rp[i];
                        baseRm[bI] = rm[i];
                        baseRs[bI] = rs[i];
                        baseQ[bI] = Q[i];
                        baseMeltHistory[bI] = meltHistory[i];
                    }
                    baseT.correctBoundaryConditions();
                }

                // Restore the bare names to base-mesh objects so write.H
                // outputs the full base mesh and writes gradT correctly.
                {
                    const fvMesh& mesh = baseMesh;
                    const volScalarField& T = baseT;

                    #include "write.H"
                }

                runTime.printExecutionTime(Info);
            }
        }
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
