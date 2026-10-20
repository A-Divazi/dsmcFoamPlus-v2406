/*---------------------------------------------------------------------------*\
  User function object: Mach number field for dsmcFoam.

    Ma = |U| / sqrt(gamma * R * T)

  where U is the mean velocity (default UMean), T the translational
  temperature (default translationalT), R = k/mass the specific gas constant.

  Usage (after dsmcFields so UMean/translationalT exist):
    dsmcMachNumber1
    {
        type            dsmcMachNumber;
        libs            (libUserDsmcFoam);
        writeControl    writeTime;
        UField          UMean;            // optional
        TField          translationalT;   // optional
        gamma           1.6666666667;     // optional (monatomic default)
        mass            66.3e-27;         // optional (argon, kg)
    }
\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_dsmcMachNumber_H
#define Foam_functionObjects_dsmcMachNumber_H

#include "fvMeshFunctionObject.H"

namespace Foam
{
namespace functionObjects
{

class dsmcMachNumber
:
    public fvMeshFunctionObject
{
    // Private Data
    word UFieldName_;
    word TFieldName_;
    scalar gamma_;
    scalar mass_;

    dsmcMachNumber(const dsmcMachNumber&) = delete;
    void operator=(const dsmcMachNumber&) = delete;

public:
    TypeName("dsmcMachNumber");

    dsmcMachNumber(const word& name, const Time& runTime, const dictionary& dict);

    virtual ~dsmcMachNumber() = default;

    virtual bool read(const dictionary&);
    virtual bool execute();
    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
