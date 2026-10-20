/*---------------------------------------------------------------------------*\
  User function object: Mach number field for dsmcFoam (implementation).
\*---------------------------------------------------------------------------*/

#include "dsmcMachNumber.H"
#include "volFields.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "OSspecific.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dsmcMachNumber, 0);

    addToRunTimeSelectionTable(functionObject, dsmcMachNumber, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dsmcMachNumber::dsmcMachNumber
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    UFieldName_("UMean"),
    TFieldName_("translationalT"),
    gamma_(5.0/3.0),
    mass_(66.3e-27)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dsmcMachNumber::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("UField", UFieldName_);
    dict.readIfPresent("TField", TFieldName_);
    dict.readIfPresent("gamma", gamma_);
    dict.readIfPresent("mass", mass_);

    return true;
}


bool Foam::functionObjects::dsmcMachNumber::execute()
{
    return true;
}


namespace Foam
{

// Load a field for the current time: from disk if present (post-process),
// otherwise from the registry (live solve).
template<class FieldType>
const FieldType* loadField
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


bool Foam::functionObjects::dsmcMachNumber::write()
{
    tmp<volVectorField> tU;
    tmp<volScalarField> tT;

    const volVectorField* UPtr = loadField(mesh_, UFieldName_, tU);
    const volScalarField* TPtr = loadField(mesh_, TFieldName_, tT);

    if (!UPtr || !TPtr)
    {
        Log << type() << ' ' << name()
            << " : missing " << UFieldName_ << " or " << TFieldName_
            << " - not calculating" << nl;

        return false;
    }

    const volVectorField& U = *UPtr;
    const volScalarField& T = *TPtr;

    const scalar R = physicoChemical::k.value()/mass_;

    volScalarField Ma
    (
        IOobject
        (
            "Ma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<scalar>("0", dimless, 0)
    );

    forAll(U.internalField(), celli)
    {
        const scalar a = Foam::sqrt(max(gamma_*R*T.internalField()[celli], SMALL));

        Ma.primitiveFieldRef()[celli] =
            (a > ROOTVSMALL) ? mag(U.internalField()[celli])/a : 0.0;
    }

    Ma.correctBoundaryConditions();

    Ma.write();

    Log << type() << ' ' << name() << " : Ma max/min : "
        << max(Ma).value() << token::SPACE << min(Ma).value() << nl;

    return true;
}


// ************************************************************************* //
