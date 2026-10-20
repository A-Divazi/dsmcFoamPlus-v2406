/*---------------------------------------------------------------------------*\
  User function object: mass flow rate through named patches (implementation).
\*---------------------------------------------------------------------------*/

#include "dsmcMassFlowRate.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include <fstream>
#include "IStringStream.H"
#include <cstdlib>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dsmcMassFlowRate, 0);

    addToRunTimeSelectionTable(functionObject, dsmcMassFlowRate, dictionary);
}
}


namespace Foam
{

//- Time directory used by the last disk read (diagnostics)
static word lastReadTime;

//- Newest written time directory of the case that contains 'name'
template<class FieldType>
fileName newestFieldPath(const fvMesh& mesh, const word& name)
{
    // local path: <case>/<time> when serial, <case>/processorN/<time> in parallel
    const fileName caseDir(mesh.time().path());
    const fileName current(mesh.time().timePath()/name);

    if (isFile(current))
    {
        return current;
    }

    // During a solve the function objects run before the fields of the
    // current time are written, so fall back to the newest time directory
    // that already contains the field (one write interval old).
    fileName newest;
    double newestValue = -1e300;

    fileNameList dirs(readDir(caseDir, fileName::DIRECTORY));

    forAll(dirs, i)
    {
        const std::string nm(dirs[i].c_str());
        char* endp = nullptr;
        const double v = std::strtod(nm.c_str(), &endp);

        // strictly numeric time directories only ("0", "0.0014", ...)
        if (endp == nm.c_str() || *endp != '\0')
        {
            continue;
        }

        const fileName f(caseDir/dirs[i]/name);

        if (isFile(f) && v > newestValue)
        {
            newestValue = v;
            newest = f;
        }
    }

    return newest;
}


//- Read exactly the given time directory (no scanning); nullptr if absent
template<class FieldType>
const FieldType* loadFieldAt
(
    const fvMesh& mesh,
    const word& name,
    const word& timeName,
    tmp<FieldType>& storage
)
{
    const fileName f(mesh.time().path()/timeName/name);

    if (!isFile(f))
    {
        return nullptr;
    }

    lastReadTime = timeName;

    storage.reset
    (
        new FieldType
        (
            IOobject
            (
                name,
                timeName,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false                   // registerObject = false
            ),
            mesh
        )
    );

    return storage.ptr();
}


template<class FieldType>
const FieldType* loadFieldDiskFirst
(
    const fvMesh& mesh,
    const word& name,
    tmp<FieldType>& storage
)
{
    const fileName path(newestFieldPath<FieldType>(mesh, name));

    if (!path.empty() && isFile(path))
    {
        // read it under the time name it was written with
        const word foundTime(path.path().name());
        lastReadTime = foundTime;

        storage.reset
        (
            new FieldType
            (
                IOobject
                (
                    name,
                    foundTime,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false                   // registerObject = false
                ),
                mesh
            )
        );

        return storage.ptr();
    }

    // Last resort: whatever is registered (post-process mode)
    return mesh.cfindObject<FieldType>(name);
}


} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::dsmcMassFlowRate::dsmcMassFlowRate
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    patchNames_(2),
    fieldName_("momentumMean"),
    rhoMFieldName_("rhoM_mix"),
    UFieldName_("UMean_mix"),
    readSource_(""),
    lastProcessed_(""),
    sections_(0),
    xTol_(5e-4)
{
    patchNames_[0] = "inlet";
    patchNames_[1] = "outlet";
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dsmcMassFlowRate::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (dict.found("patches"))
    {
        patchNames_ = dict.lookup("patches");
    }
    dict.readIfPresent("field", fieldName_);
    dict.readIfPresent("rhoMField", rhoMFieldName_);
    dict.readIfPresent("UField", UFieldName_);
    dict.readIfPresent("sections", sections_);
    dict.readIfPresent("sectionsTol", xTol_);

    if (!sections_.empty())
    {
        buildSections();
    }

    return true;
}


void Foam::functionObjects::dsmcMassFlowRate::buildSections()
{
    secFaces_.clear();
    secSf_.clear();
    secOwner_.clear();
    secNames_.clear();

    const fvMesh& mesh = mesh_;
    const label nInternal = mesh.nInternalFaces();
    const labelUList& own = mesh.owner();
    const vectorField Sf(mesh.faceAreas());
    const vectorField C(mesh.faceCentres());

    for (const scalar x0 : sections_)
    {
        DynamicList<label> fl;
        DynamicList<vector> sf;
        DynamicList<label> oc;

        for (label f = 0; f < nInternal; ++f)
        {
            if (mag(C[f].x() - x0) < xTol_)
            {
                fl.append(f);
                sf.append(Sf[f]);
                oc.append(own[f]);
            }
        }

        secFaces_.append(fl);
        secSf_.append(sf);
        secOwner_.append(oc);
        std::ostringstream os;
        os << "x=" << x0;
        secNames_.append(os.str());
    }
}


bool Foam::functionObjects::dsmcMassFlowRate::execute()
{
    // dsmcFoam+ writes its fields after the function objects of the same write
    // step have run, so wait until a new time directory with the fields exists
    // and process that. Each write time is processed exactly once.
    const fileName path(newestFieldPath<volScalarField>(mesh_, rhoMFieldName_));

    if (path.empty() || !isFile(path))
    {
        return true;
    }

    const word foundTime(path.path().name());

    if (foundTime == lastProcessed_)
    {
        return true;
    }

    if (process(foundTime))
    {
        lastProcessed_ = foundTime;
    }

    return true;
}


bool Foam::functionObjects::dsmcMassFlowRate::process(const word& foundTime)
{
    tmp<volVectorField> tMom;
    const volVectorField* momPtr = nullptr;

    // Read exactly the time directory being processed; if it is not complete
    // yet (the directory can exist before its fields are written) report it as
    // unprocessed so the next execute() call retries it.
    {
        tmp<volScalarField> tRho0;
        tmp<volVectorField> tU0;
        const volScalarField* r0 =
            loadFieldAt(mesh_, rhoMFieldName_, foundTime, tRho0);
        const volVectorField* u0 =
            loadFieldAt(mesh_, UFieldName_, foundTime, tU0);

        if (!r0 || !u0)
        {
            return false;
        }
    }

    // dsmcFoam+ has no momentum field: build it as rhoM*U from the fields it
    // does write (density and mean velocity).
    tmp<volScalarField> tRho;
    tmp<volVectorField> tU;
    if (!momPtr)
    {
        const volScalarField* rhoPtr =
            loadFieldAt(mesh_, rhoMFieldName_, foundTime, tRho);
        const volVectorField* uPtr =
            loadFieldAt(mesh_, UFieldName_, foundTime, tU);

        if (rhoPtr && uPtr)
        {
            tMom.reset
            (
                new volVectorField
                (
                    IOobject
                    (
                        fieldName_,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    (*rhoPtr)*(*uPtr)
                )
            );

            momPtr = tMom.ptr();
        }

        readSource_ = lastReadTime;
    }

    if (!momPtr)
    {
        Log << type() << ' ' << name() << " : no " << fieldName_ << " (nor "
            << rhoMFieldName_ << " and " << UFieldName_
            << ") field found - not calculating" << nl;

        return false;
    }

    const volVectorField& mom = *momPtr;


    scalar total = 0;
    scalarList mdots;
    wordList labels;

    if (!sections_.empty())
    {
        // Measure the flux through interior cross-section planes (established
        // flow) — avoids the inflow boundary-cell deficit that biases the
        // inlet patch measurement.
        labels = secNames_;
        mdots.setSize(secFaces_.size(), 0.0);

        forAll(secFaces_, s)
        {
            const labelList& fl = secFaces_[s];
            const List<vector>& sf = secSf_[s];
            const labelList& oc = secOwner_[s];

            scalar mdot = 0;
            forAll(fl, f)
            {
                mdot += mom.internalField()[oc[f]] & sf[f];
            }

            reduce(mdot, sumOp<scalar>());

            mdots[s] = mdot;
        }

        // For interior sections both fluxes are positive (flow +x); the
        // "net" mass-conservation residual is their difference (0 at steady).
        total = (mdots.size() >= 2) ? (mdots.last() - mdots[0]) : 0.0;
    }
    else
    {
        // Measure at the boundary patches (inlet/outlet).
        labelList patchIDs(patchNames_.size(), -1);
        forAll(patchNames_, i)
        {
            patchIDs[i] = mesh_.boundaryMesh().findPatchID(patchNames_[i]);
        }

        labels = patchNames_;
        mdots.setSize(patchIDs.size(), 0.0);

        forAll(patchIDs, i)
        {
            const label pi = patchIDs[i];

            if (pi < 0)
            {
                Log << type() << ' ' << name() << " : patch " << patchNames_[i]
                    << " not found" << nl;
                continue;
            }

            const polyPatch& patch = mesh_.boundaryMesh()[pi];
            const labelList& faceCells = patch.faceCells();
            const vectorField& Sf = patch.faceAreas();

            scalar mdot = 0;
            forAll(faceCells, facei)
            {
                mdot += mom.internalField()[faceCells[facei]] & Sf[facei];
            }

            reduce(mdot, sumOp<scalar>());

            mdots[i] = mdot;
            total += mdot;
        }
    }

    // Append to a single CSV (master processor only).
    // globalPath() = case root (path() would give .../processor0 in parallel).
    if (Pstream::master())
    {
        const fileName dir(mesh_.time().globalPath()/"postProcessing");
        mkDir(dir);
        const fileName csv(dir/"massFlow.csv");

        const bool first = !isFile(csv);

        std::ofstream os(csv.c_str(), std::ios::app);
        if (first)
        {
            os << "#time";
            for (const word& lb : labels) os << ' ' << lb;
            os << " net" << std::endl;
        }
        os << Foam::readScalar(IStringStream(foundTime.c_str())());
        for (const scalar m : mdots) os << ' ' << m;
        os << ' ' << total << std::endl;
    }

    Log << type() << ' ' << name() << " : ";
    forAll(labels, i)
    {
        Log << labels[i] << '=' << mdots[i] << " kg/s ";
    }
    Log << " net=" << total << " kg/s" << nl;

    return true;
}


bool Foam::functionObjects::dsmcMassFlowRate::write()
{
    return true;
}


// ************************************************************************* //
