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
    boundaryMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "boundaryMeasurements.H"
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
boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{}


//- Construct from mesh, cloud and boolean (dsmcFoam)
boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const bool& dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    typeIds_(cloud_.typeIdList().size(),-1),
    rhoNIntBF_(),
    rhoNElecBF_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    mccSpeciesBF_(),
    momentumBF_(),
    UMeanBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    vibrationalEBF_(),
    electronicEBF_(),
    qBF_(),
    fDBF_()
{
    
    forAll(typeIds_, i)
    {
        typeIds_[i] = i;
    }
    
    rhoNIntBF_.setSize(typeIds_.size());
    rhoNElecBF_.setSize(typeIds_.size());
    rhoNBF_.setSize(typeIds_.size());
    rhoMBF_.setSize(typeIds_.size());
    linearKEBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());
    momentumBF_.setSize(typeIds_.size());
    UMeanBF_.setSize(typeIds_.size());
    rotationalEBF_.setSize(typeIds_.size());
    rotationalDofBF_.setSize(typeIds_.size());
    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    qBF_.setSize(typeIds_.size());
    fDBF_.setSize(typeIds_.size());
    
    
//     forAll(rhoNIntBF_, i)
//     {
//         const polyPatch& patch = mesh_.boundaryMesh()[i];
//         rhoNIntBF_[i].setSize(patch.size(),0.0);
//         rhoNElecBF_[i].setSize(patch.size(),0.0);
//     }
    
    forAll(rhoNBF_, i)
    {        
        rhoNIntBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoNElecBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoMBF_[i].setSize(mesh_.boundaryMesh().size());
        linearKEBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        momentumBF_[i].setSize(mesh_.boundaryMesh().size());
        UMeanBF_[i].setSize(mesh_.boundaryMesh().size());
        rotationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        rotationalDofBF_[i].setSize(mesh_.boundaryMesh().size());
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        qBF_[i].setSize(mesh_.boundaryMesh().size());
        fDBF_[i].setSize(mesh_.boundaryMesh().size());
        //****//
        
        forAll(rhoNBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            rhoNIntBF_[i][j].setSize(patch.size(),0.0);
            rhoNElecBF_[i][j].setSize(patch.size(),0.0);
            rhoNBF_[i][j].setSize(patch.size(),0.0);
            rhoMBF_[i][j].setSize(patch.size(),0.0);
            linearKEBF_[i][j].setSize(patch.size(),0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(),0.0);
            momentumBF_[i][j].setSize(patch.size(),vector::zero);
            UMeanBF_[i][j].setSize(patch.size(),vector::zero);
            rotationalEBF_[i][j].setSize(patch.size(),0.0);
            rotationalDofBF_[i][j].setSize(patch.size(),0.0);
            vibrationalEBF_[i][j].setSize(patch.size(),0.0);
            electronicEBF_[i][j].setSize(patch.size(),0.0);
            qBF_[i][j].setSize(patch.size(),0.0);
            fDBF_[i][j].setSize(patch.size(),vector::zero);
            //****//
        }
    }
    
    
    
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundaryMeasurements::~boundaryMeasurements()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void boundaryMeasurements::clean()
{
    //- clean geometric fields
    
    forAll(rhoNBF_, i)
    {
        forAll(rhoNBF_[i], j)
        {
            rhoNBF_[i][j] = 0.0;
            rhoMBF_[i][j] = 0.0;
            linearKEBF_[i][j] = 0.0;
            mccSpeciesBF_[i][j] = 0.0;
            momentumBF_[i][j] = vector::zero;
            UMeanBF_[i][j] = vector::zero;
            rotationalEBF_[i][j] = 0.0;
            rotationalDofBF_[i][j] = 0.0;
            vibrationalEBF_[i][j] = 0.0;
            electronicEBF_[i][j] = 0.0;
            qBF_[i][j] = 0.0;
            fDBF_[i][j] = vector::zero;
            //****//
            rhoNIntBF_[i][j] = 0.0;
            rhoNElecBF_[i][j] = 0.0;
        }
    }
}


void boundaryMeasurements::updateFields
(
    dsmcParcel& p
)
{}

void boundaryMeasurements::autoMap
(
    const mapPolyMesh& mapper
)
{
    // store old BFs
    // Outer list is typeId, middle is the patch, inner is the local patch face
//    List<List<scalarField>> rhoNIntBFOld_ = rhoNIntBF_;
//    List<List<scalarField>> rhoNElecBFOld_ = rhoNElecBF_;
//    List<List<scalarField>> rhoNBFOld_ = rhoNBF_;
//    List<List<scalarField>> rhoMBFOld_ = rhoMBF_;
//    List<List<scalarField>> linearKEBFOld_ = linearKEBF_;
//    List<List<scalarField>> mccSpeciesBFOld_ = mccSpeciesBF_;
//    List<List<vectorField>> momentumBFOld_ = momentumBF_;
//    List<List<vectorField>> UMeanBFOld_ = UMeanBF_;
//    List<List<scalarField>> rotationalEBFOld_ = rotationalEBF_;
//    List<List<scalarField>> rotationalDofBFOld_ = rotationalDofBF_;
//    List<List<scalarField>> vibrationalEBFOld_ = vibrationalEBF_;
//    List<List<scalarField>> electronicEBFOld_ = electronicEBF_;
//    List<List<scalarField>> qBFOld_ = qBF_;
//    List<List<vectorField>> fDBFOld_ = fDBF_;

    forAll(typeIds_, i)
    {
        typeIds_[i] = i;
    }

    rhoNIntBF_.setSize(typeIds_.size());
    rhoNElecBF_.setSize(typeIds_.size());
    rhoNBF_.setSize(typeIds_.size());
    rhoMBF_.setSize(typeIds_.size());
    linearKEBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());
    momentumBF_.setSize(typeIds_.size());
    UMeanBF_.setSize(typeIds_.size());
    rotationalEBF_.setSize(typeIds_.size());
    rotationalDofBF_.setSize(typeIds_.size());
    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    qBF_.setSize(typeIds_.size());
    fDBF_.setSize(typeIds_.size());

    forAll(rhoNBF_, i) {
        rhoNIntBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoNElecBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoMBF_[i].setSize(mesh_.boundaryMesh().size());
        linearKEBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        momentumBF_[i].setSize(mesh_.boundaryMesh().size());
        UMeanBF_[i].setSize(mesh_.boundaryMesh().size());
        rotationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        rotationalDofBF_[i].setSize(mesh_.boundaryMesh().size());
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        qBF_[i].setSize(mesh_.boundaryMesh().size());
        fDBF_[i].setSize(mesh_.boundaryMesh().size());
        //****//

        forAll(rhoNBF_[i], j) {
            const polyPatch &patch = mesh_.boundaryMesh()[j];
            rhoNIntBF_[i][j].setSize(patch.size(), 0.0);
            rhoNElecBF_[i][j].setSize(patch.size(), 0.0);
            rhoNBF_[i][j].setSize(patch.size(), 0.0);
            rhoMBF_[i][j].setSize(patch.size(), 0.0);
            linearKEBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
            momentumBF_[i][j].setSize(patch.size(), vector::zero);
            UMeanBF_[i][j].setSize(patch.size(), vector::zero);
            rotationalEBF_[i][j].setSize(patch.size(), 0.0);
            rotationalDofBF_[i][j].setSize(patch.size(), 0.0);
            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            qBF_[i][j].setSize(patch.size(), 0.0);
            fDBF_[i][j].setSize(patch.size(), vector::zero);
            //****//
            }
        }
    }

//    const labelList& faceMap = mapper.faceMap();
//    const labelList& reverseFaceMap = mapper.reverseFaceMap();
//
//    // Storage for any master faces. These will be the original faces
//    // on the coarse cell that get split into four (or rather the
//    // master face gets modified and three faces get added from the master)
//    // Estimate number of faces created
//
//    bitSet masterFaces(nFaces());
//
//    forAll(faceMap, facei)
//    {
//        const label oldFacei = faceMap[facei];
//
//        if (oldFacei >= 0)
//        {
//            const label masterFacei = reverseFaceMap[oldFacei];
//
//            if (masterFacei < 0)
//            {
//                FatalErrorInFunction
//                << "Problem: should not have removed faces"
//                << " when refining."
//                << nl << "face:" << facei << endl
//                << abort(FatalError);
//            }
//            else if (masterFacei != facei)
//            {
//                masterFaces.set(masterFacei);
//            }
//        }
//    }
//
//}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
