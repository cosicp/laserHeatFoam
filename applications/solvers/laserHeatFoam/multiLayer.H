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

    // Helper: collect face-neighbours of `frontier` that are not yet in
    // `excludeA` or `excludeB`.
    const auto haloStep =
        [&cellCells]
        (
            const labelHashSet& frontier,
            const labelHashSet& excludeA,
            const labelHashSet& excludeB
        ) -> labelHashSet
        {
            labelHashSet next(frontier.size());

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

            if (next.empty())
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

            if (next.empty())
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
    Info<< "    layer " << i << " : " << layerCells[i].size() << " cells"
        << endl;
    totalLayerCells += layerCells[i].size();
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

Info<< "    total active cells across all layers: " << totalLayerCells
    << " / " << mesh.nCells() << " base mesh cells";

if (totalLayerCells < mesh.nCells())
{
    Info<< " (" << (mesh.nCells() - totalLayerCells)
        << " cells will remain inactive throughout the run)";
}

Info<< nl << endl;


// ************************************************************************* //
