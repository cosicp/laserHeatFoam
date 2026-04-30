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

    Provides the following identifiers in enclosing scope:
        nLayers           label
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
const scalar layerDuration = multiLayerDict.get<scalar>("layerDuration");
const scalar layerInitialT = multiLayerDict.get<scalar>("layerInitialT");
const word exposedFacesPatch =
    multiLayerDict.lookupOrDefault<word>("exposedFacesPatch", word::null);

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

    labelHashSet layer0(2*bp.size());
    forAll(bpFaceCells, i)
    {
        layer0.insert(bpFaceCells[i]);
    }
    layerCells[0] = std::move(layer0);

    labelHashSet covered(layerCells[0]);
    const labelListList& cellCells = mesh.cellCells();

    for (label i = 1; i < nLayers; ++i)
    {
        labelHashSet next(layerCells[i-1].size());

        for (const label cellI : layerCells[i-1])
        {
            const labelList& nbrs = cellCells[cellI];

            forAll(nbrs, n)
            {
                const label nbr = nbrs[n];

                if (!covered.found(nbr))
                {
                    next.insert(nbr);
                }
            }
        }

        if (next.empty())
        {
            FatalErrorInFunction
                << "Layer " << i << " has no cells. The mesh appears to"
                << " have fewer than " << nLayers << " layers reachable"
                << " from patch '" << baseplatePatch << "'."
                << exit(FatalError);
        }

        covered |= next;
        layerCells[i] = std::move(next);
    }
}

label totalLayerCells = 0;
Info<< "Multi-layer: per-layer cell counts (face-neighbour growth from '"
    << baseplatePatch << "'):" << endl;

forAll(layerCells, i)
{
    Info<< "    layer " << i << " : " << layerCells[i].size() << " cells"
        << endl;
    totalLayerCells += layerCells[i].size();
}

Info<< "    total active cells across all layers: " << totalLayerCells
    << " / " << mesh.nCells() << " base mesh cells";

if (totalLayerCells < mesh.nCells())
{
    Info<< " (" << (mesh.nCells() - totalLayerCells)
        << " cells will remain inactive throughout the run)";
}

Info<< nl << endl;


// ************************************************************************* //
