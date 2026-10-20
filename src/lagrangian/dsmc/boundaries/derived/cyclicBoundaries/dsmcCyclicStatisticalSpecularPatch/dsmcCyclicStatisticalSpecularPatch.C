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

Description

\*---------------------------------------------------------------------------*/

#include "dsmcCyclicStatisticalSpecularPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(dsmcCyclicStatisticalSpecularPatch, 1);

    addToRunTimeSelectionTable(dsmcCyclicBoundary, dsmcCyclicStatisticalSpecularPatch, dictionary);


    void dsmcCyclicStatisticalSpecularPatch::readProperties()
    {
        Tr0_ = (readScalar(propsDict_.lookup("Tr")));
        velocity_ = (propsDict_.get<vector>("velocity"));
        omegaLim_ = (acos(1 - Foam::sqrt(Tr0_)));
    }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
    dsmcCyclicStatisticalSpecularPatch::dsmcCyclicStatisticalSpecularPatch
            (
                    Time& t,
                    const polyMesh& mesh,
                    dsmcCloud& cloud,
                    const dictionary& dict
            )
            :
            dsmcCyclicBoundary(t, mesh, cloud, dict),
            propsDict_(dict.subDict(typeName + "Properties")),
            Tr0_(readScalar(propsDict_.lookup("Tr"))),
            omegaLim_(acos( 1 - Foam::sqrt(Tr0_))),
            velocity_(vector(propsDict_.lookup("velocity"))),
            nReflections_(0),
            nRejections_(0),
            meanOmega_(0.0)
    {
        writeInTimeDir_ = false;
        writeInCase_ = false;
    }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    dsmcCyclicStatisticalSpecularPatch::~dsmcCyclicStatisticalSpecularPatch()
    {}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
    scalar dsmcCyclicStatisticalSpecularPatch::calculateTransmission
            (
                    const scalar& omega
            )
    {
        scalar xTr1 = ((1 - sqrt(Tr0_)) / 2.) * (1. / cos(omega) -1);
//    Info << "xTr1 = " << xTr1 << endl;
        scalar Tr = ((1 - ((1 - sqrt(Tr0_)) + 2.0 * xTr1))) * sqrt(Tr0_);
//    Info << "Tr = " << Tr << endl;
        return Tr;
    }

    void dsmcCyclicStatisticalSpecularPatch::calculateProperties()
    {
        label nReflections = nReflections_;
        label nRejections = nRejections_;
        scalar meanOmega = meanOmega_;

        if(Pstream::parRun())
        {
            reduce(nReflections, sumOp<label>());
            reduce(nRejections, sumOp<label>());
            reduce(meanOmega, sumOp<scalar>());
        }
        if (nReflections > 0 || nRejections > 0) {
            Info << "Number of Reflections = " << nReflections
                 << ", Number of Rejections = " << nRejections
                 << " Ratio Relfections/(Reflections+Rejections) = "
                 << scalar(nReflections) / scalar(nReflections + nRejections)
                 << ", Mean Omega = " << (meanOmega *180.0 / pi) / scalar(nReflections + nRejections)
                 << endl;
        }
    }

    void dsmcCyclicStatisticalSpecularPatch::initialConfiguration()
    {}



    void dsmcCyclicStatisticalSpecularPatch::controlMol
            (
                    dsmcParcel& mol,
                    dsmcParcel::trackingData& td
            )
    {
        const label& faceI = mol.face();
        vector nF = mesh_.faceAreas()[faceI];
        nF /= mag(nF);

        scalar cosOmega = (nF & mol.U())/
                          (mag(nF) * mag(mol.U())) + VSMALL; // radian angle between nF and mol.U()
        scalar cosOmega2 = cosOmega * cosOmega;
        scalar omega = acos(cosOmega);

        meanOmega_ += omega;

//    Info << "omega = " << omega * 180 / pi << " in Limits [" << omegaLim_*180/ pi << ", " << (2.0*pi -omegaLim_)*180 / pi << "?" << endl;
        if(omega > omegaLim_)
        {
            p_ =  scalar(0.0);
        }
        else
        {
            p_ = calculateTransmission(omega);
        }

//    Info<< "parcel to reflect at pos: "
//        << mol.position() << ", nF: " << nF
//        << " old velocity = " << mol.U()
//        << " faceI = " << faceI
//        << " omega = " << omega * 180 / pi
//        << " transmission prob.= " << p_
//        << endl;

        Random& rndGen = cloud_.rndGen();
        scalar pRandom = rndGen.sample01<double>();

//    Info << "pRandom = " << pRandom << "; " << "p_ = " << p_ << endl;

        if (pRandom >= p_) // reflect
        {
            scalar Un = mol.U() & nF;
            mol.U() -= 2.0*Un*nF;
            mol.U() += velocity_;

            nReflections_++;
            reflected_ = true;


//         Pout<< "Reflected!!: mol at pos: "
//             << mol.position() << ", nF: " << nF
//             << " new velocity: " << mol.U()
//             << endl;
        }
        else
        {
            nRejections_++;
            reflected_ = false;
        }

    }



    void dsmcCyclicStatisticalSpecularPatch::output
            (
                    const fileName& fixedPathName,
                    const fileName& timePath
            )
    {

    }


    void dsmcCyclicStatisticalSpecularPatch::updateProperties(const dictionary& newDict)
    {
        //- the main properties should be updated first
        updateBoundaryProperties(newDict);

        readProperties();
    }


    void dsmcCyclicStatisticalSpecularPatch::autoBoundaryModelMap
            (
                    const mapPolyMesh& mapper
            )
    {
    }

    bool dsmcCyclicStatisticalSpecularPatch::reflected()
    {
        return reflected_;
    }


} // End namespace Foam

// ************************************************************************* //
