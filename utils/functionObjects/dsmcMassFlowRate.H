/*---------------------------------------------------------------------------*\
  User function object: mass flow rate through named patches, for dsmcFoam.

    mdot = integral of (momentum . n) dA  over each patch, with
    momentum = rho*U (default: momentumMean, the time-averaged field).

  Appends one line per writeTime to  postProcessing/massFlow.csv :
      time  mdot_<patch1>  mdot_<patch2> ...  net

  Sign: outward face normals, so inflow (x=0) is negative, outflow positive.

  Usage (after fieldAverage so momentumMean exists):
    dsmcMassFlowRate1
    {
        type            dsmcMassFlowRate;
        libs            (libUserDsmcFoam);
        writeControl    writeTime;
        patches         (inlet outlet);
        field           momentumMean;
    }
\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_dsmcMassFlowRate_H
#define Foam_functionObjects_dsmcMassFlowRate_H

#include "fvMeshFunctionObject.H"

namespace Foam
{
namespace functionObjects
{

class dsmcMassFlowRate
:
    public fvMeshFunctionObject
{
    // Private Data
    wordList patchNames_;
    word fieldName_;
    //- Fields used to build momentum (=rhoM*U) when the requested vector
    //  field is not available (e.g. dsmcFoam+, which writes rhoM_<x>, UMean_<x>)
    word rhoMFieldName_;
    word UFieldName_;
    //- time directory the fields were actually read from (diagnostics)
    mutable word readSource_;
    //- newest time directory already processed (polling in execute())
    word lastProcessed_;
    scalarList sections_;      // empty = use patches; else interior x-planes
    scalar xTol_;
    // Precomputed section-face data (for each section: faces -> Sf -> owner cell)
    List<List<label>> secFaces_;
    List<List<vector>> secSf_;
    List<labelList> secOwner_;
    List<word> secNames_;

    //- Build internal face-plane data for the requested sections
    void buildSections();

    dsmcMassFlowRate(const dsmcMassFlowRate&) = delete;
    void operator=(const dsmcMassFlowRate&) = delete;

public:
    TypeName("dsmcMassFlowRate");

    dsmcMassFlowRate(const word& name, const Time& runTime, const dictionary& dict);

    virtual ~dsmcMassFlowRate() = default;

    virtual bool read(const dictionary&);
    virtual bool execute();

    bool process(const word& foundTime);

    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
