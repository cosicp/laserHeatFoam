/*---------------------------------------------------------------------------*\
Description
    Refresh refineIndicator on the base mesh.

    meltHistory is binary (0 in unconsolidated cells, 1 in cells that
    have ever crossed the liquidus). dynamicRefineFvMesh refines on
    cell-centre values inside (lowerRefineLevel, upperRefineLevel), so
    a raw 0/1 field never falls in any (0.1, 0.9) band.

    We build a smoothed surrogate:

        ms = fvc::average-smoothed(meltHistory),  nSmoothing sweeps
        refineIndicator = 4 * ms * (1 - ms)

    The 4*x*(1-x) map peaks at 1 in the interface band (ms ~ 0.5) and
    falls to 0 deep in either consolidated (ms -> 1) or unconsolidated
    (ms -> 0) regions, so a single one-sided unrefineLevel correctly
    coarsens BOTH ends of the binary jump while refining the interface.

    fvc::average evaluates a face-area-weighted neighbour-cell average
    via surfaceInterpolation, which goes through the field boundary
    conditions including processor patches — so each smoothing sweep is
    parallel-correct without any explicit halo handling here.

    nSmoothing is captured at startup from constant/dynamicMeshDict
    (key "indicatorSmoothing", default 3) and lives as a local in the
    solver — keeping the setting co-located with the refinement
    parameters in the dict, even though it's applied here on the base
    mesh and not by dynamicRefineFvMesh itself.
\*---------------------------------------------------------------------------*/

{
    volScalarField ms(baseMeltHistory);

    for (label sweep = 0; sweep < indicatorSmoothing; ++sweep)
    {
        ms = fvc::average(ms);
        ms.correctBoundaryConditions();
    }

    refineIndicator = scalar(4) * ms * (scalar(1) - ms);
}


// ************************************************************************* //
