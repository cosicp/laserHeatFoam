    if (runTime.writeTime())
    {
    
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
