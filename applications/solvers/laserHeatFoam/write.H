if (runTime.writeTime())
{
    // Derived field written only at output times to keep time-step cost low.
    volVectorField gradT
    (
        IOobject
        (
            "gradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(T)
    );

    runTime.write();
}
