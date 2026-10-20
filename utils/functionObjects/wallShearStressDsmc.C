/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    User function object: wall shear stress for dsmcFoam (implementation).

    See wallShearStressDsmc.H for usage.
\*---------------------------------------------------------------------------*/

#include "wallShearStressDsmc.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallShearStressDsmc, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        wallShearStressDsmc,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallShearStressDsmc::wallShearStressDsmc
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fDFieldName_("fDMean"),
    tauFieldName_("wallShearStress"),
    tauMagFieldName_("wallShearStressMag")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallShearStressDsmc::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("fDField", fDFieldName_);
    dict.readIfPresent("tauField", tauFieldName_);
    dict.readIfPresent("tauMagField", tauMagFieldName_);

    return true;
}


bool Foam::functionObjects::wallShearStressDsmc::execute()
{
    return true;
}


bool Foam::functionObjects::wallShearStressDsmc::write()
{
    tmp<volVectorField> tfD;
    const volVectorField* fDPtr = nullptr;

    // Prefer reading the field fresh from disk (post-process mode, where the
    // registry can hold a single time's field). During a live solve the field
    // for the current time is not yet on disk, so fall back to the registry.
    const fileName fDFile(mesh_.time().timePath()/fDFieldName_);

    if (isFile(fDFile))
    {
        tfD.reset
        (
            new volVectorField
            (
                IOobject
                (
                    fDFieldName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            )
        );

        fDPtr = tfD.ptr();
    }
    else
    {
        fDPtr = mesh_.cfindObject<volVectorField>(fDFieldName_);
    }

    if (!fDPtr)
    {
        Log << type() << ' ' << name() << " : no " << fDFieldName_
            << " field found - not calculating" << nl;

        return false;
    }

    const volVectorField& fD = *fDPtr;

    // Volume fields so the results appear in ParaView/Tecplot (surface
    // fields are ignored by the standard OpenFOAM readers). The internal
    // field is zero; only the wall boundary values are meaningful.
    volVectorField tau
    (
        IOobject
        (
            tauFieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<vector>(fD.dimensions(), Zero),
        "calculated"
    );

    volScalarField tauMag
    (
        IOobject
        (
            tauMagFieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<scalar>(fD.dimensions(), Zero),
        "calculated"
    );

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];

        if (isA<wallPolyPatch>(patch))
        {
            const vectorField n(patch.faceAreas()/mag(patch.faceAreas()));

            const vectorField& f = fD.boundaryField()[patchi];

            const vectorField t(f - (f & n)*n);   // tangential projection

            tau.boundaryFieldRef()[patchi] = t;

            tauMag.boundaryFieldRef()[patchi] = mag(t);

            // Also map the shear into the wall-adjacent CELLS, so it is
            // visible and extractable (line/slice through the near-wall
            // cells or surface colouring) in ParaView/Tecplot.
            for (label facei = 0; facei < t.size(); ++facei)
            {
                const label celli = mesh_.owner()[patch.start() + facei];

                tau.primitiveFieldRef()[celli] = t[facei];

                tauMag.primitiveFieldRef()[celli] = mag(t[facei]);
            }
        }
    }

    tau.write();
    tauMag.write();

    Log << type() << ' ' << name() << " : " << tauFieldName_
        << " mag max/min : "
        << max(mag(tau.boundaryField())) << token::SPACE
        << min(mag(tau.boundaryField())) << nl;

    return true;
}


// ************************************************************************* //
