/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    cellMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "cellMeasurements.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from mesh and cloud 
cellMeasurements::cellMeasurements
(
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{}


//- Construct from mesh, cloud and boolean (dsmcFoam)
cellMeasurements::cellMeasurements
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const bool& dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    collisionSeparation_(),
    nColls_()
{
    collisionSeparation_.setSize(mesh.nCells());
    nColls_.setSize(mesh.nCells());   
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cellMeasurements::~cellMeasurements()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cellMeasurements::clean()
{
    //- clean geometric fields
    
    collisionSeparation_ = 0.0;
    nColls_ = 0.0;
}


void cellMeasurements::updateFields
(
    dsmcParcel& p
)
{
}

void cellMeasurements::autoMap
(
    const mapPolyMesh& mapper
)
{
    const labelList &cellMap = mapper.cellMap();
    const labelList &reverseCellMap = mapper.reverseCellMap();

    // Each split cell becomes original + 7 additional
    labelList nSubCells(mapper.nOldCells(), 0);

    forAll(cellMap, celli)
    {
        const label oldCelli = cellMap[celli];
        if (oldCelli >= 0 && reverseCellMap[oldCelli] >= 0)
        {
            // Found master cell.
            nSubCells[oldCelli]++;
        }
    }

    scalarField correctedCollisionSeparation_;
    scalarField correctedNColls_;

    correctedCollisionSeparation_.setSize(mesh_.nCells());
    correctedNColls_.setSize(mesh_.nCells());

    // Correct any split cells
    forAll(cellMap, celli)
    {
        const label oldCelli = cellMap[celli];
        const scalar fRatio = mesh_.V()[celli] / mesh_.V0()[oldCelli];
        Pout << "Cell ratio V/V0 = " << fRatio << endl;

        if (oldCelli >= 0 && nSubCells[oldCelli] == 8) {
            correctedCollisionSeparation_[celli] = collisionSeparation()[oldCelli]; // intensive
            correctedNColls_[celli] = nColls()[oldCelli] * fRatio; // extensive
        }

        const auto &cellsFromCells = mapper.cellsFromCellsMap();
        for (const auto &s: cellsFromCells) {
            // Reset unsplit cell
            const label celli = s.index();
            // sum up extensive vars
            const labelList &oldCells = s.masterObjects();
            for (const label oldCelli: oldCells) {
                correctedCollisionSeparation_[celli] = collisionSeparation()[oldCelli]; // intensive
                correctedNColls_[celli] += nColls()[oldCelli]; // extensive
            }
        }
    }

    collisionSeparation() = correctedCollisionSeparation_;
    nColls() = correctedNColls_;

}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
