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
    Heat conduction solver for Laser Powder Bed Fusion (LPBF) with a volumetric heat source.

    \heading Solver details
    The solver solves the heat conduction equation for a scalar quantity, T.  The
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
#include "cmath"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote(
        "Heat conduction solver for Laser Powder Bed Fusion (LPBF) with a volumetric heat source.");

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating temperature distribution\n"
         << endl;

        while (runTime.loop())
        {
            Info << "Time = " << runTime.timeName() << nl << endl;

            // Store previous liquid fraction (for post-processing / optional source-form)
            fLold = fL;

            #include "updateThermo.H"
            #include "updateLaser.H"
            #include "DiffusionNo.H"

            // Apparent heat capacity latent heat model (Proell et al. 2020 / 2023):
            // - liquid fraction fL ramps linearly between Ts and Tl
            // - cpEff = cp + L * dfL/dT, where dfL/dT = 1/(Tl-Ts) in mushy zone else 0
            // Notes:
            // - L is specific latent heat [J/kg]
            // - rho*cpEff enters the transient term (volumetric heat capacity)

            const dimensionedScalar dTmush
            (
                "dTmush", dimTemperature, max(SMALL, (Tl - Ts).value())
            );

            fL = min
            (
                max((T - Ts)/dTmush, dimensionedScalar("zero", dimless, 0)),
                dimensionedScalar("one",  dimless, 1)
            );

            const volScalarField mushyMask(pos(T - Ts)*pos(Tl - T));
            cpEff = cp + mushyMask*(L/dTmush);

            // Consolidated fraction (Proell 2020, eq. 23):
            // rc = max fL ever reached — irreversible powder->melt transition.
            rc = max(rc, fL);

            // Phase fractions (Proell 2020, eqs. 24-26):
            //   rp = 1 - rc   (powder fraction)
            //   rm = fL        (melt fraction)
            //   rs = rc - fL   (solid fraction)

            // Phase fractions (Proell 2020, eqs. 24-26):
            rp = scalar(1) - rc;
            rm = fL;
            rs = rc - fL;

            // Phase-dependent thermal conductivity (Proell 2020, eq. 27):
            //   k = rp*kp + rm*km + rs*ks
            k = rp*kp + rm*km + rs*ks;

            // Binary melt marker (1 if cell ever fully melted)
            meltHistory = max(meltHistory, pos(T - Tl));

            volScalarField rhoCpEff
            (
                IOobject
                (
                    "rhoCpEff",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rho*cpEff
            );

            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix TEqn(
                    fvm::ddt(rhoCpEff, T)
                    - fvm::laplacian(k, T)
                    == Q
                );

                // fvOptions.constrain(TEqn);
                TEqn.solve();
                // fvOptions.correct(T);
            }

            Info << "T after sol: " << max(T) << endl;
            #include "write.H"
            runTime.printExecutionTime(Info);
        }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
