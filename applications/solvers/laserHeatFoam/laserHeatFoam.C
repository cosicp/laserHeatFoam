/*---------------------------------------------------------------------------*\
  =========       E          |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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
    Heat conduction solver for Laser Powder Bed Fusion (LPBF) with a volumetric
    heat source.

    \heading Solver details
    The solver solves the heat conduction equation for a scalar quantity, T. The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( \alpha \grad T \right) + Q
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        alpha | Thermal diffusivity
        Q     | Volumetric heat source
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "interpolationTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Heat conduction solver for Laser Powder Bed Fusion (LPBF) "
        "with a volumetric heat source."
    );

#include "postProcess.H"
#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"

    simpleControl simple(mesh);

#include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Preserve the previous liquid fraction before thermo updates.
        fLold = fL;

        #include "updateThermo.H"
        #include "updateLaser.H"
        #include "DiffusionNo.H"
        #include "updateCpeff.H"
        #include "updateFractions.H"

        // Rebuilding the effective conductivity from the current phase mix.
        k = rp*kp + rm*km + rs*ks;

        // Solve the transient heat equation with an apparent heat capacity.
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                rho*cpEff*fvm::ddt(T) - fvm::laplacian(k, T) == Q
            );

            TEqn.solve();
        }

        Info<< "T after sol (max, full domain): " << max(T) << endl;

        // Optional diagnostic below the laser surface, useful for monitoring
        // substrate heating without post-processing the full field.
        // if (substrateProbeDepth > SMALL)
        // {
        //     scalar maxTsubstrate = -GREAT;

        //     forAll(T, cellI)
        //     {
        //         if (mesh.C()[cellI].z() <= laserSurfaceZ - substrateProbeDepth)
        //         {
        //             maxTsubstrate = max(maxTsubstrate, T[cellI]);
        //         }
        //     }

        //     reduce(maxTsubstrate, maxOp<scalar>());

        //     Info<< "T after sol (substrate, z<=laserZ-"
        //         << substrateProbeDepth << "): " << maxTsubstrate
        //         << " K" << endl;
        // }

        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
