/*---------------------------------------------------------------------------*\
Description
    Rebuild the active-cell hash set from the persistent layerID field.

    A cell is active iff 0 <= layerID[cellI] <= currentLayer. Called once
    per layer activation (so cells with layerID == currentLayer come on
    line) and after every base-mesh refinement event (so the hash uses
    the post-refinement cell labels rather than stale pre-refinement
    ones).

    Tolerant to floating-point representation: layerID is stored as a
    scalar but only ever takes integer (or -1) values. The lround +
    bounds check accepts e.g. 5.0 and 4.999999 alike.
\*---------------------------------------------------------------------------*/

inline void rebuildActiveCells
(
    const volScalarField& layerID,
    const Foam::label currentLayer,
    Foam::labelHashSet& activeCells
)
{
    activeCells.clear();
    activeCells.resize(2*layerID.size());

    const Foam::scalarField& fld = layerID.primitiveField();
    forAll(fld, cellI)
    {
        const Foam::label k = Foam::label(std::lround(fld[cellI]));

        if (k >= 0 && k <= currentLayer)
        {
            activeCells.insert(cellI);
        }
    }
}


// ************************************************************************* //
