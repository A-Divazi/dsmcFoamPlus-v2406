/*---------------------------------------------------------------------------*\
  User function object: density-gradient fields for dsmcFoam.

    gradRho    = nabla(rho)         (volVectorField)
    gradRhoMag = |nabla(rho)|       (volScalarField, "numerical Schlieren")

  rho defaults to rhoMMean (time-averaged mass density from fieldAverage).

  Usage (after fieldAverage so rhoMMean exists):
    dsmcGradRho1
    {
        type            dsmcGradRho;
        libs            (libUserDsmcFoam);
        writeControl    writeTime;
        rhoField        rhoMMean;     // optional
        gradField       gradRho;      // optional
        magField        gradRhoMag;   // optional
    }
\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_dsmcGradRho_H
#define Foam_functionObjects_dsmcGradRho_H

#include "fvMeshFunctionObject.H"

namespace Foam
{
namespace functionObjects
{

class dsmcGradRho
:
    public fvMeshFunctionObject
{
    // Private Data
    word rhoFieldName_;
    word gradFieldName_;
    word magFieldName_;

    dsmcGradRho(const dsmcGradRho&) = delete;
    void operator=(const dsmcGradRho&) = delete;

public:
    TypeName("dsmcGradRho");

    dsmcGradRho(const word& name, const Time& runTime, const dictionary& dict);

    virtual ~dsmcGradRho() = default;

    virtual bool read(const dictionary&);
    virtual bool execute();
    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
