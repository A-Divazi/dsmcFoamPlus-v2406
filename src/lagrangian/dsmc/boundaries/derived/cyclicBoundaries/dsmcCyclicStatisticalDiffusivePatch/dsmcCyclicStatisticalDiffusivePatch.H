/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    dsmcCyclicStatisticalDiffusivePatch

Description

    Statistical description for mesh structures, with given (mesh) parameters gridWidthB, gridLengthK.
    Models the reflection probability as a function of the angle beween face normal vector and velocity
    vector of the particle.


SourceFiles
    dsmcCyclicStatisticalDiffusivePatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcCyclicStatisticalDiffusivePatch_H
#define dsmcCyclicStatisticalDiffusivePatch_H

#include "dsmcCyclicBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcCyclicStatisticalDiffusivePatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcCyclicStatisticalDiffusivePatch
:
	public dsmcCyclicBoundary
{

private:

    // Private data
    dictionary propsDict_;


    // Private Member Functions

    void readProperties();

    bool reflected_;

    scalar p_;

    scalar Tr0_;
    scalar omegaLim_;

    scalar temperature_;
    vector velocity_;

    label nReflections_;
    label nRejections_;


public:

    //- Runtime type information
    TypeName("dsmcCyclicStatisticalDiffusivePatch");

    // Constructors

        //- Construct from components
        dsmcCyclicStatisticalDiffusivePatch
		(
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~dsmcCyclicStatisticalDiffusivePatch();


    // Member Functions
        scalar calculateTransmission
        (
            const scalar& omega
        );


        void calculateProperties();

        //- create an initial configuration
        void initialConfiguration();

        //- control the molecules
        void controlMol
        (
            dsmcParcel& mol,
            dsmcParcel::trackingData& td
        );

        bool reflected();

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);


        // Access

        //- map
        void autoBoundaryModelMap(const mapPolyMesh& mapper);





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
