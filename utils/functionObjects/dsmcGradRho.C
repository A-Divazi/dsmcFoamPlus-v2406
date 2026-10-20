/*---------------------------------------------------------------------------*\
  User function object: density-gradient fields for dsmcFoam (implementation).
\*---------------------------------------------------------------------------*/

#include "dsmcGradRho.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "dimensionSets.H"
#include "addToRunTimeSelectionTable.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dsmcGradRho, 0);

    addToRunTimeSelectionTable(functionObject, dsmcGradRho, dictionary);
}
}


namespace Foam
{

// Load a field for the current time: from disk if present (post-process),
// otherwise from the registry (live solve).
template<class FieldType>
const FieldType* loadFieldDiskFirst
(
    const fvMesh& mesh,
    const word& name,
    tmp<FieldType>& storage
)
{
    const fileName path(mesh.time().timePath()/name);

    if (isFile(path))
    {
        storage.reset
        (
            new FieldType
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );

        return storage.ptr();
    }

    return mesh.cfindObject<FieldType>(name);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dsmcGradRho::dsmcGradRho
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    rhoFieldName_("rhoMMean"),
    gradFieldName_("gradRho"),
    magFieldName_("gradRhoMag")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dsmcGradRho::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("rhoField", rhoFieldName_);
    dict.readIfPresent("gradField", gradFieldName_);
    dict.readIfPresent("magField", magFieldName_);

    return true;
}


bool Foam::functionObjects::dsmcGradRho::execute()
{
    return true;
}


bool Foam::functionObjects::dsmcGradRho::write()
{
    tmp<volScalarField> tRho;
    const volScalarField* rhoPtr =
        loadFieldDiskFirst(mesh_, rhoFieldName_, tRho);

    if (!rhoPtr)
    {
        Log << type() << ' ' << name() << " : no " << rhoFieldName_
            << " field found - not calculating" << nl;

        return false;
    }

    const volScalarField& rho = *rhoPtr;

    volVectorField gradRho
    (
        IOobject
        (
            gradFieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<vector>("0", rho.dimensions()/dimLength, Zero)
    );

    gradRho = fvc::grad(rho);

    volScalarField gradRhoMag
    (
        IOobject
        (
            magFieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<scalar>("0", rho.dimensions()/dimLength, Zero)
    );

    gradRhoMag = mag(gradRho);

    gradRho.write();
    gradRhoMag.write();

    Log << type() << ' ' << name() << " : " << gradFieldName_
        << " mag max/min : "
        << max(mag(gradRho)).value() << token::SPACE
        << min(mag(gradRho)).value() << nl;

    return true;
}


// ************************************************************************* //
