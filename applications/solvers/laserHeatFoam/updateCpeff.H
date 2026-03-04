forAll(mesh.cells(), cellI)
{
    scalar Tc = T[cellI];

    if (Tc > Ts && Tc < Tl)
    {
        // Info << "Cell " << cellI << " is in mushy zone. Adjusting cpEff.\n" << endl;
        cpEff[cellI] = cp[cellI] + L/(Tl - Ts);
    }
    else
    {
        // Info << "Cell " << cellI << " is outside mushy zone. cpEff = cp.\n" << endl;
        cpEff[cellI] = cp[cellI];
    }

}