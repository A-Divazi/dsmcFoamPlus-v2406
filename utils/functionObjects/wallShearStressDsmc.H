/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    User function object: wall shear stress for dsmcFoam.

    Derives the wall shear stress vector (and its magnitude) from the wall
    boundary values of a DSMC force-density field (default: fDMean):

        tau_wall = fD - (fD . n) n

    where n is the wall face normal. For DSMC, fD boundary values are the
    particle momentum-transfer rate per unit area (normal part = pressure,
    tangential part = shear).

    Usage in system/controlDict.functions (after fieldAverage):

        wallShearStress1
        {
            type            wallShearStressDsmc;
            libs            (libUserDsmcFoam);
            writeControl    writeTime;
            fDField         fDMean;             // optional, default fDMean
            tauField        wallShearStress;    // optional output names
            tauMagField     wallShearStressMag;
        }
\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_wallShearStressDsmc_H
#define Foam_functionObjects_wallShearStressDsmc_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class wallShearStressDsmc Declaration
\*---------------------------------------------------------------------------*/

class wallShearStressDsmc
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of the DSMC force-density field (default: fDMean)
        word fDFieldName_;

        //- Output field names
        word tauFieldName_;
        word tauMagFieldName_;


    // Private Member Functions

        //- No copy construct
        wallShearStressDsmc(const wallShearStressDsmc&) = delete;

        //- No copy assignment
        void operator=(const wallShearStressDsmc&) = delete;


public:

    //- Runtime type information
    TypeName("wallShearStressDsmc");


    // Constructors

        //- Construct from name, Time and dictionary
        wallShearStressDsmc
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~wallShearStressDsmc() = default;


    // Member Functions

        //- Read the function-object dictionary
        virtual bool read(const dictionary&);

        //- Execute the function-object operations, currently does nothing
        virtual bool execute();

        //- Write the function-object results
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
