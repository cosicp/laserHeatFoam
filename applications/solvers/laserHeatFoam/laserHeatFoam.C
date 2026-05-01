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
#include "dynamicFvMesh.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "interpolationTable.H"
#include "mathematicalConstants.H"
#include "fvMeshSubset.H"
#include "rebuildActiveCells.H"

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
    #include "createDynamicFvMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "multiLayer.H"

    // Number of Jacobi smoothing sweeps applied to meltHistory when
    // building refineIndicator. Co-located with the rest of the AMR
    // configuration in constant/dynamicMeshDict; default 3 sweeps if
    // absent. Read once at startup; the dict itself is owned by the
    // dynamicFvMesh.
    const label indicatorSmoothing =
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        ).getOrDefault<label>("indicatorSmoothing", 3);

    // Re-bind base-mesh fields under explicit names so we can shadow the
    // bare names with sub-mesh equivalents inside the layer loop without
    // losing access to the base versions.
    dynamicFvMesh& baseMesh = mesh;
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
    volScalarField& baseLayerID = layerID;

    // fvMeshSubset's subMesh internals reference the base mesh's
    // current topology; after a base-mesh refinement event the subset
    // has to be rebuilt from scratch (fresh fvMeshSubset), not just
    // reset(). Hold by autoPtr so we can swap it after each AMR step.
    autoPtr<fvMeshSubset> subMesherPtr(new fvMeshSubset(baseMesh));
    labelHashSet activeCells(2*baseMesh.nCells());

    // dynamicRefineFvMesh::mapFields reads the base-mesh old-cell-volume
    // field V0 when remapping conserved fields across a refinement
    // event. fvMesh only stores V0 inside updateMesh() if VPtr_ has been
    // realised (fvMesh.C:977 "if (VPtr_) storeOldVol(...)"). We solve T
    // on the sub-mesh and never touch base-mesh V in the steady state,
    // so without this nudge VPtr_ stays null and the first refinement
    // event aborts with "V0 is not available". Materialising V() once
    // at startup is sufficient — fvMesh keeps it alive and storeOldVol
    // takes over from there.
    (void)baseMesh.V();

    Info<< nl
        << "Simulating " << nLayers << " layer(s); each layer runs for "
        << layerDuration << " s." << nl
        << "Total simulated time = " << nLayers*layerDuration << " s." << nl
        << "(controlDict endTime is overridden per-layer.)" << nl << endl;

    Info<< "Calculating temperature distribution" << nl << endl;

    for (label layerI = 0; layerI < nLayers; ++layerI)
    {
        // Initialise temperature on cells being activated this iteration
        // (those with layerID == layerI). Already-active cells keep their
        // running thermal history. Using baseLayerID as the source rather
        // than a stored labelHashSet means the activation step is robust
        // to any AMR that occurred during previous layers — children of
        // a layer-K cell inherit layerID == K, so they are still picked
        // up here by their post-refinement labels.
        const scalarField& lidFld = baseLayerID.primitiveField();
        label nNewCells = 0;
        forAll(lidFld, cellI)
        {
            if (label(std::lround(lidFld[cellI])) == layerI)
            {
                baseT[cellI] = layerInitialT;
                ++nNewCells;
            }
        }
        baseT.correctBoundaryConditions();

        Info<< "Activating layer " << layerI
            << " (" << nNewCells << " new cells)" << endl;

        // Grow the active set to include this layer and rebuild the
        // sub-mesh.
        rebuildActiveCells(baseLayerID, layerI, activeCells);
        subMesherPtr->reset(activeCells, exposedPatchID);

        // Extend simulation end time to cover this layer.
        runTime.setEndTime((layerI + 1)*layerDuration);

        Info<< "    sub-mesh: " << subMesherPtr->subMesh().nCells()
            << " cells; layer end time: " << runTime.endTime().value()
            << " s" << nl << endl;

        // Time loop for this layer. The sub-mesh and its derived fields
        // are constructed inside the loop body so they are rebuilt
        // freshly when AMR has changed the base-mesh topology.
        while (runTime.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;

            // dynamicRefineFvMesh::update() reads refineIndicator from
            // the base mesh and refines/unrefines accordingly. All
            // registered volScalarFields (T, rho, cp, rc, meltHistory,
            // layerID, refineIndicator, ...) are mapped automatically.
            const bool meshChanged = baseMesh.update();

            if (meshChanged)
            {
                Info<< "    AMR: base mesh now has "
                    << baseMesh.nCells() << " cells" << endl;

                // Cell labels have been remapped; rebuild activeCells
                // from the (mapped) layerID and rebuild the subset
                // helper from the refined base mesh.
                rebuildActiveCells(baseLayerID, layerI, activeCells);
                subMesherPtr.reset(new fvMeshSubset(baseMesh));
                subMesherPtr->reset(activeCells, exposedPatchID);
            }

            fvMeshSubset& subMesher = subMesherPtr();

            // Sub-mesh scope: the bare names mesh, T, rho, ... are
            // shadowed here so that the existing update*.H, DiffusionNo.H
            // and write.H includes operate on the sub-mesh without
            // modification.
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
            // write.H produces complete output for the whole domain
            // and so the next mesh.update() sees an up-to-date
            // baseMeltHistory when computing refineIndicator.
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
            } // end sub-mesh scope (sub-mesh fields destroyed here)

            // Refresh the AMR driver field from the (just-updated)
            // baseMeltHistory. This must be done after the sub-mesh
            // scope closes so its sub-mesh "meltHistory" object is no
            // longer in the registry shadowing the base-mesh field.
            #include "updateRefineIndicator.H"

            // Write the full-base-mesh output. write.H references
            // mesh and T; both already resolve to the base versions
            // (the sub-mesh shadows are out of scope).
            #include "write.H"

            runTime.printExecutionTime(Info);
        }
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
