/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dsmcParcel.H"
#include "dsmcCloud.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dsmcParcel::move
(
    dsmcCloud& cloud,
    dsmcParcel::trackingData& td,
    const scalar trackTime
)
{
    dsmcParcel& p =
            static_cast<dsmcParcel&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    if(stuckToWall() == 0)
    {
        const polyMesh& mesh = cloud.pMesh();
        const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

        if(newParcel() == 1)
        {
            Random& rndGen(cloud.rndGen());
            //stepFraction() = rndGen.scalar01();
            stepFraction() = rndGen.sample01<double>();
            newParcel() = 0;
        }
        
        
        //scalar tEnd = (1.0 - stepFraction())*trackTime;
        //const scalar dtMax = tEnd;

                    
        // For reduced-D cases, the velocity used to track needs to be
        // constrained, but the actual U_ of the parcel must not be
        // altered or used, as it is altered by patch interactions an
        // needs to retain its 3D value for collision purposes.
        vector Utracking = U_;
            
        while (td.keepParticle && !td.switchProcessor && p.stepFraction() < 1)
        {           
            Utracking = U_;
            vector d = vector::zero;

            if(!cloud.axisymmetric())
            {
                // Apply correction to position for reduced-D cases, 
                // but not axisymmetric cases
                //meshTools::constrainToMeshCentre(mesh, position());
                
                // Apply correction to velocity to constrain tracking for
                // reduced-D cases,  but not axisymmetric cases
                meshTools::constrainDirection
                (
                    mesh, mesh.solutionD(), Utracking
                );
                // Deviation from the mesh centre for reduced-D cases
                vector d = p.deviationFromMeshCentre();

            }

            const scalar f = 1 - p.stepFraction();

            // Set the Lagrangian time-step
            //scalar dt = min(dtMax, tEnd);

            //dt *= trackToFace(position() + dt*Utracking-d, td, true);
            p.trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);

        }
    }

    return td.keepParticle;
}


bool Foam::dsmcParcel::hitPatch
(
    dsmcCloud&,
    dsmcParcel::trackingData&
)
{
    return false;
}

void Foam::dsmcParcel::hitProcessorPatch
(
    dsmcCloud&,
    dsmcParcel::trackingData& td
)
{
    td.switchProcessor = true;
}

void Foam::dsmcParcel::hitWallPatch
(
    dsmcCloud& cloud,
    dsmcParcel::trackingData& td
)
{
    const label wppIndex = this->patch();
    //-find which patch has been hit

    const label& patchModelId = cloud.boundaries().
    patchToModelIds()[wppIndex];

    // apply a boundary model when a molecule collides with this poly patch
    cloud.boundaries().
    patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
}

void Foam::dsmcParcel::hitWedgePatch
(
        dsmcCloud& cloud,
        dsmcParcel::trackingData& td
)
{
    hitSymmetryPatch(cloud, td);
}


void Foam::dsmcParcel::hitCyclicPatch
(
        dsmcCloud& cloud,
        dsmcParcel::trackingData& td
)
{
    const label wppIndex = this->patch();
    //-find which patch has been hit

    const label& cyclicBModelId = cloud.boundaries()
            .cyclicBoundaryToModelIds()[wppIndex];


    // apply a boundary model when a molecule collides with this poly patch
    if(cyclicBModelId!=-1)
    {
        cloud.boundaries().
                cyclicBoundaryModels()[cyclicBModelId]->controlMol(*this, td);
        if (!cloud.boundaries().
                cyclicBoundaryModels()[cyclicBModelId]->reflected())
        {
            particle::hitCyclicPatch(cloud, td);
        }

    }
    else // i.e., no cyclicModel defined
    {
        particle::hitCyclicPatch(cloud, td);
    }
}


void Foam::dsmcParcel::hitCyclicAMIPatch
(
        dsmcCloud& cloud,
        dsmcParcel::trackingData& td,
        const vector& displacement
)
{
    const label wppIndex = this->patch();
    //-find which patch has been hit

    const label& cyclicBModelId = cloud.boundaries()
            .cyclicBoundaryToModelIds()[wppIndex];


    // apply a boundary model when a molecule collides with this poly patch
    if(cyclicBModelId!=-1)
    {
        cloud.boundaries().
                cyclicBoundaryModels()[cyclicBModelId]->controlMol(*this, td);
        if (!cloud.boundaries().
                cyclicBoundaryModels()[cyclicBModelId]->reflected())
        {
            particle::hitCyclicAMIPatch(cloud, td, displacement);
        }

    }
    else // i.e., no cyclicModel defined
    {
        particle::hitCyclicAMIPatch(cloud, td, displacement);
    }
}


void Foam::dsmcParcel::hitCyclicACMIPatch
        (
                dsmcCloud& cloud,
                dsmcParcel::trackingData& td,
                const vector& displacement
        )
{
    const label wppIndex = this->patch();
    //-find which patch has been hit

    const label& cyclicBModelId = cloud.boundaries()
            .cyclicBoundaryToModelIds()[wppIndex];


    // apply a boundary model when a molecule collides with this poly patch
    if(cyclicBModelId!=-1)
    {
        cloud.boundaries().
                cyclicBoundaryModels()[cyclicBModelId]->controlMol(*this, td);
        if (!cloud.boundaries().
                cyclicBoundaryModels()[cyclicBModelId]->reflected())
        {
            particle::hitCyclicACMIPatch(cloud, td, displacement);
        }

    }
    else // i.e., no cyclicModel defined
    {
        particle::hitCyclicACMIPatch(cloud, td, displacement);
    }
}

 //template<class ParcelType>
void Foam::dsmcParcel::transformProperties
(
    const tensor& T
)
{
   particle::transformProperties(T);
   U_ = transform(T, U_);
}


 //template<class ParcelType>
void Foam::dsmcParcel::transformProperties
(
    const vector& separation
)
{
  particle::transformProperties(separation);
}

//void Foam::dsmcParcel::correctTet()
//{
//    // See note in correctAfterParallelTransfer for tetPti addressing ...
//    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;
//}
//
//void Foam::dsmcParcel::setTopology
//(
//    label cellI,
//    label faceI
//)
//{
//    celli_ = cellI;
//    facei_ = faceI;
//}
//
//void Foam::dsmcParcel::setCoordinates
//(
//    scalar a,
//    scalar b,
//    scalar c,
//    scalar d
//)
//{
//    coordinates_.a() = a;
//    coordinates_.b() = b;
//    coordinates_.c() = c;
//    coordinates_.d() = d;
//}

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "dsmcParcelIO.C"

// ************************************************************************* //
