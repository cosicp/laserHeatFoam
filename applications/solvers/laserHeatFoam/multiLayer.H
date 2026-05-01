/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Description
    Multi-layer activation control for laserHeatFoam.

    Reads constant/multiLayerProperties and pre-computes which base-mesh
    cells belong to each layer, growing outward from a user-named
    baseplate patch using face-neighbour topology. Layer 0 is the set of
    cells with a face on the baseplate patch; layer k (k > 0) is the set
    of face-neighbours of layer k - 1 that have not yet been visited.

    Each user-level layer can correspond to M >= 1 mesh cells in the
    thickness direction (cellsPerLayer); the layer is then built as M
    successive face-neighbour growth steps starting from the baseplate
    patch (for layer 0) or from the boundary of the previous layer.

    Provides the following identifiers in enclosing scope:
        nLayers           label
        cellsPerLayer     label   (mesh cells through the thickness, >= 1)
        layerDuration     scalar  [s]
        layerInitialT     scalar  [K]
        exposedPatchID    label   (-1 -> default oldInternalFaces patch)
        layerCells        List<labelHashSet>, size == nLayers

\*---------------------------------------------------------------------------*/

IOdictionary multiLayerDict
(
    IOobject
    (
        "multiLayerProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const word baseplatePatch = multiLayerDict.get<word>("baseplatePatch");
const label nLayers = multiLayerDict.get<label>("nLayers");
const label cellsPerLayer =
    multiLayerDict.lookupOrDefault<label>("cellsPerLayer", 1);
const scalar layerDuration = multiLayerDict.get<scalar>("layerDuration");
const scalar layerInitialT = multiLayerDict.get<scalar>("layerInitialT");
const word exposedFacesPatch =
    multiLayerDict.lookupOrDefault<word>("exposedFacesPatch", word::null);

if (cellsPerLayer < 1)
{
    FatalErrorInFunction
        << "cellsPerLayer must be >= 1, got " << cellsPerLayer
        << exit(FatalError);
}

label exposedPatchID = -1;
if (!exposedFacesPatch.empty())
{
    exposedPatchID = mesh.boundaryMesh().findPatchID(exposedFacesPatch);

    if (exposedPatchID < 0)
    {
        FatalErrorInFunction
            << "exposedFacesPatch '" << exposedFacesPatch
            << "' not found in base mesh boundary." << nl
            << "Available patches: " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    Info<< "Multi-layer: exposed sub-mesh faces will be added to existing"
        << " base patch '" << exposedFacesPatch << "' (id "
        << exposedPatchID << ")." << endl;
}
else
{
    Info<< "Multi-layer: no exposedFacesPatch specified; exposed faces"
        << " will be placed on the synthetic 'oldInternalFaces' patch"
        << " (zeroGradient-equivalent default)." << endl;
}

const label baseplatePatchID =
    mesh.boundaryMesh().findPatchID(baseplatePatch);

if (baseplatePatchID < 0)
{
    FatalErrorInFunction
        << "baseplatePatch '" << baseplatePatch
        << "' not found in base mesh boundary." << nl
        << "Available patches: " << mesh.boundaryMesh().names()
        << exit(FatalError);
}

// Per-layer cell sets, built by face-neighbour growth from the baseplate.
List<labelHashSet> layerCells(nLayers);

{
    const polyPatch& bp = mesh.boundaryMesh()[baseplatePatchID];
    const labelUList& bpFaceCells = bp.faceCells();
    const labelListList& cellCells = mesh.cellCells();

    // Cells already assigned to some layer.
    labelHashSet covered(2*bp.size()*cellsPerLayer);

    // Seed for the next layer: cells from which the next halo step grows.
    // Initially this is the set of cells with a face on the baseplate.
    labelHashSet seed(2*bp.size());
    forAll(bpFaceCells, i)
    {
        seed.insert(bpFaceCells[i]);
    }

    // Reusable scalar marker used by haloStep to broadcast the frontier
    // across processor patches. Allocating it once outside the BFS loop
    // avoids the per-step IOobject + field construction cost.
    volScalarField frontierMarker
    (
        IOobject
        (
            "frontierMarker",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Helper: collect face-neighbours of `frontier` that are not yet in
    // `excludeA` or `excludeB`. Two stages so the BFS is correct in
    // parallel:
    //   1. local expansion via mesh.cellCells() (intra-rank face graph).
    //   2. cross-processor expansion: paint the frontier into a
    //      volScalarField, evaluate boundary conditions so processor
    //      patches exchange their owner-cell values, then any local
    //      cell whose coupled-patch face has neighbour value > 0.5 is
    //      a face-neighbour of a frontier cell on the other rank.
    // For serial runs the second stage is a no-op (no coupled patches
    // carry remote owners) so the original behaviour is preserved.
    const auto haloStep =
        [&mesh, &cellCells, &frontierMarker]
        (
            const labelHashSet& frontier,
            const labelHashSet& excludeA,
            const labelHashSet& excludeB
        ) -> labelHashSet
        {
            labelHashSet next(frontier.size());

            // 1. Local face-neighbour expansion.
            for (const label cellI : frontier)
            {
                const labelList& nbrs = cellCells[cellI];

                forAll(nbrs, n)
                {
                    const label nbr = nbrs[n];

                    if (!excludeA.found(nbr) && !excludeB.found(nbr))
                    {
                        next.insert(nbr);
                    }
                }
            }

            // 2. Cross-processor expansion (no-op in serial).
            if (Pstream::parRun())
            {
                scalarField& mFld = frontierMarker.primitiveFieldRef();
                mFld = scalar(0);
                for (const label cellI : frontier)
                {
                    mFld[cellI] = scalar(1);
                }
                frontierMarker.correctBoundaryConditions();

                forAll(mesh.boundary(), patchI)
                {
                    if (!mesh.boundary()[patchI].coupled())
                    {
                        continue;
                    }

                    const fvPatchScalarField& pf =
                        frontierMarker.boundaryField()[patchI];
                    const labelUList& fcs =
                        mesh.boundary()[patchI].faceCells();

                    forAll(fcs, i)
                    {
                        if (pf[i] > scalar(0.5))
                        {
                            const label cellI = fcs[i];
                            if
                            (
                                !excludeA.found(cellI)
                             && !excludeB.found(cellI)
                            )
                            {
                                next.insert(cellI);
                            }
                        }
                    }
                }
            }

            return next;
        };

    for (label k = 0; k < nLayers; ++k)
    {
        // For k > 0, replace the seed with the first halo step from the
        // previous layer's cells (so that the seed lies just above the
        // last fully-assigned layer).
        if (k > 0)
        {
            labelHashSet next =
                haloStep(layerCells[k-1], covered, labelHashSet());

            // Only fatal if no rank found a new cell. Per-rank emptiness
            // is normal: a layer may live entirely on a different
            // processor.
            if (returnReduce(next.size(), sumOp<label>()) == 0)
            {
                FatalErrorInFunction
                    << "Layer " << k << " has no cells. The mesh has"
                    << " fewer than " << nLayers << "*" << cellsPerLayer
                    << " cells in the thickness direction reachable"
                    << " from patch '" << baseplatePatch << "'."
                    << exit(FatalError);
            }

            seed = std::move(next);
        }

        // Build this layer: seed + (cellsPerLayer - 1) more halo steps.
        labelHashSet thisLayer(seed);
        labelHashSet frontier(seed);

        for (label step = 1; step < cellsPerLayer; ++step)
        {
            labelHashSet next = haloStep(frontier, covered, thisLayer);

            // Globally fatal only: per-rank emptiness is normal.
            if (returnReduce(next.size(), sumOp<label>()) == 0)
            {
                FatalErrorInFunction
                    << "Layer " << k << " halo step " << step
                    << " of " << (cellsPerLayer - 1) << " produced no"
                    << " new cells. The mesh has fewer than "
                    << cellsPerLayer << " cells in the thickness"
                    << " direction at this layer."
                    << exit(FatalError);
            }

            thisLayer |= next;
            frontier = std::move(next);
        }

        covered |= thisLayer;
        layerCells[k] = std::move(thisLayer);
    }
}

label totalLayerCells = 0;
Info<< "Multi-layer: per-layer cell counts (face-neighbour growth from '"
    << baseplatePatch << "', " << cellsPerLayer
    << " cell(s) per layer through the thickness):" << endl;

forAll(layerCells, i)
{
    // In parallel each rank holds only its slice; sum for the global
    // count that the user actually wants to see.
    const label sz =
        returnReduce(layerCells[i].size(), sumOp<label>());
    Info<< "    layer " << i << " : " << sz << " cells" << endl;
    totalLayerCells += sz;
}

// Populate the persistent base-mesh layerID field from the BFS sets.
// This is the only place layerCells[] is consumed; once written, layerID
// is mapped automatically across base-mesh refinement events and is the
// single source of truth for layer membership thereafter.
forAll(layerCells, i)
{
    for (const label cellI : layerCells[i])
    {
        layerID[cellI] = scalar(i);
    }
}
layerID.correctBoundaryConditions();

const label globalNCells = returnReduce(mesh.nCells(), sumOp<label>());

Info<< "    total active cells across all layers: " << totalLayerCells
    << " / " << globalNCells << " base mesh cells";

if (totalLayerCells < globalNCells)
{
    Info<< " (" << (globalNCells - totalLayerCells)
        << " cells will remain inactive throughout the run)";
}

Info<< nl << endl;


// ************************************************************************* //
