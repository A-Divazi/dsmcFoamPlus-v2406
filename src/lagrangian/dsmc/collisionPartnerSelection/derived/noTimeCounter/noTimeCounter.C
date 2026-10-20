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
    noTimeCounter

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounter.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable
(collisionPartnerSelection, noTimeCounter, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noTimeCounter::noTimeCounter
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    infoCounter_(0)
//     propsDict_(dict.subDict(typeName + "Properties"))
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounter::~noTimeCounter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounter::initialConfiguration()
{

}

void noTimeCounter::collide()
{
    //Pout<< "This is noTimeCounter::collide()" << nl << endl;
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    List<DynamicList<label> > subCells(8);

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;

    const List<DynamicList<dsmcParcel*> >& cellOccupancy =
                                                cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();

    //Pout<< "Search collission partners " << nl << endl;

    forAll(cellOccupancy, cellI)
    {
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);
        const scalar& cellVolume = mesh.cellVolumes()[cellI];

        label nC(cellParcels.size());

        if (nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            point cC = mesh.cellCentres()[cellI];
            //Pout << "for all cellParcells.." << endl;
            forAll(cellParcels, i)
            {
                const dsmcParcel& p = *cellParcels[i];

                //Pout << "calc rel. position of particle " << i << " of cell " << cellI << "..." << endl;
                //cloud_.dumpParticlePositions();

                vector relPos = p.position() - cC;

                //Pout << "rel. position: " << p.position() << endl;

                label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                //Pout << "subCell=" << subCell << endl;

                //Pout << "subCells[suCell]=" << subCells[subCell] << "and append " << i << endl;
                subCells[subCell].append(i);

                whichSubCell[i] = subCell;
                //Pout << "==" << nl << endl;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //Pout << "Selecting Pairs in Cell " << cellI << endl;
            scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];
            if (sigmaTcRMax < VSMALL)
            {
                Pout << "sigmaTcRMax < VSMALL: " << sigmaTcRMax << nl << endl;
            }
            
            scalar selectedPairs = 0.0;

            if (cloud_.collisionSelectionRemainder()[cellI]<VSMALL) // prevent numerical issues with very small cSRs
            {
                cloud_.collisionSelectionRemainder()[cellI] = 0;
                Pout << "cSR < VSMALL" << nl << endl;
            }
            else if (cloud_.collisionSelectionRemainder()[cellI]>1)
            {
                cloud_.collisionSelectionRemainder()[cellI] = 1-VSMALL;
                Pout << "cSR > 1" << nl << endl;
            }
               
            if(cloud_.axisymmetric())
            {
                //Pout<< "for axis symmetric simulation." << nl << endl;
                scalar RWF = 0.0;
                scalar nMols = 0.0;
                
                forAll(cellParcels, i)
                {
                    const dsmcParcel& p = *cellParcels[i];
                    
                    scalar radius = sqrt(sqr(p.position().y()) 
                                        + sqr(p.position().z()));

                    RWF += 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                    
                    nMols += 1.0;
                }

                RWF /= nMols;
                
//                 const point& cC = mesh.cellCentres()[cellI];
                
//                 scalar radius = cC.y();
                
//                 RWF = 1.0 + cloud_.maxRWF()
//                                         *(radius/cloud_.radialExtent());
                selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*cloud_.nParticle()*RWF*sigmaTcRMax*deltaT
                /cellVolume;
            }
            else
            {
                selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*cloud_.nParticle()*sigmaTcRMax*deltaT
                /cellVolume;
            }

//            Pout    <<"cSR=" << cloud_.collisionSelectionRemainder()[cellI] <<  "nC=" << nC << " FN=" << cloud_.nParticle() << " sigmaTcRMax="
//                    << sigmaTcRMax << " deltaT=" << deltaT << " Vcell="<< cellVolume << endl;
//            Pout << "Calc = " << cloud_.collisionSelectionRemainder()[cellI]
//                                 + 0.5*nC*(nC - 1)*cloud_.nParticle()*sigmaTcRMax*deltaT
//                                   /cellVolume << endl;

            //Pout << "selectedPairs = " << selectedPairs << endl;
            label nCandidates(selectedPairs);
            //Pout << "nCandidates(selectedPairs) " << nCandidates << endl;

            cloud_.collisionSelectionRemainder()[cellI] = selectedPairs 
                                                                - nCandidates;

            //Pout << "collisionSelectionRemainder()[cellI] " << cloud_.collisionSelectionRemainder()[cellI] << endl;
            collisionCandidates += nCandidates;
            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                //label candidateP = rndGen_.integer(0, nC - 1);
                label candidateP = 0 + lrand48() % (nC - 0);
                //Pout << "candidateP " << candidateP << nl << endl;

                // Declare the second collision candidate
                label candidateQ = -1;

                const List<label>& subCellPs = 
                                        subCells[whichSubCell[candidateP]];

                label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.
                    do
                    {
                        //candidateQ = subCellPs[rndGen_.integer(0, nSC - 1)];
                        candidateQ = 0 + (lrand48() % (nSC - 0));


                    } while (candidateP == candidateQ);
                }
                else
                {

                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        //candidateQ = rndGen_.integer(0, nC - 1);
                        candidateQ = 0 + (lrand48() % (nC - 0));

                    } while (candidateP == candidateQ);
                }
                //Pout << "Candidate selected" << nl << endl;
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = rndGen_.integer(0, nC-1);

                // // Select a possible second collision candidate
                // label candidateQ = rndGen_.integer(0, nC-1);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = rndGen_.integer(0, nC-1);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                dsmcParcel& parcelP = *cellParcels[candidateP];
                dsmcParcel& parcelQ = *cellParcels[candidateQ];

                label chargeP = -2;
                label chargeQ = -2;

                chargeP = cloud_.constProps(parcelP.typeId()).charge();
                chargeQ = cloud_.constProps(parcelQ.typeId()).charge();
                
                //do not allow electron-electron collisions
                
                if(!(chargeP == -1 && chargeQ == -1))
                {

                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP,
                        parcelQ
                    );
                    

                    // Update the maximum value of sigmaTcR stored, but use the
                    // initial value in the acceptance-rejection criteria 
                    // because the number of collision candidates selected was 
                    // based on this


                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    //if ((sigmaTcR/sigmaTcRMax) > rndGen_.scalar01())
                    if ((sigmaTcR/sigmaTcRMax) > rndGen_.sample01<double>())
                    {
                        // chemical reactions

                        // find which reaction model parcel p and q should use
                        label rMId = cloud_.reactions().returnModelId
                                                        (parcelP, parcelQ);

                        if(rMId != -1)
                        {
                                //Pout << "chemical reaction:" << endl;
                                cloud_.reactions().reactions()[rMId]->reaction
                                (
                                    parcelP,
                                    parcelQ
                                );                                   
                            // if reaction unsuccessful use conventional
                            //    collision model
                            if(cloud_.reactions().reactions()[rMId]->relax())
                            {
                                //Pout << "neutral collision: reaction unsuccessful." << endl;
                                cloud_.binaryCollision().collide
                                (
                                    parcelP,
                                    parcelQ,
                                    cellI
                                );
                            }
                        }
                        // if reaction model not found, use conventional 
                        // collision model
                        else 
                        {
                            //Pout << "neutral collision: no reaction model found." << endl;
                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }
                        //Pout << "collision successful!" << nl << endl;

                        collisions++;
                    }
                }
            }

            //Pout << "collisions done." << nl << endl;
        }
    }

    //Pout << "Done." << nl << endl;
    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    //Pout << "Correct Boundary condition.." << endl;

    cloud_.sigmaTcRMax().correctBoundaryConditions();
    
    infoCounter_++;
        
    if(infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisionCandidates)
        {
            Info<< "    Collisions                      = "
                << collisions << nl
                << "    Acceptance rate                 = "
                << scalar(collisions)/scalar(collisionCandidates) << nl
                << endl;
                
            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No collisions" << endl;
            
            infoCounter_ = 0;
        }
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
