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

#include "dsmcCyclicForwardScatterStatisticalDiffusivePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "../grid.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
const std::vector<double> x_vals = TransConstants::Tr0;
const std::vector<double> y_vals = TransConstants::omega;
const std::vector<std::vector < double>> z_vals = TransConstants::funcTr;
namespace Foam
{

defineTypeNameAndDebug(dsmcCyclicForwardScatterStatisticalDiffusivePatch, 1);

addToRunTimeSelectionTable(dsmcCyclicBoundary, dsmcCyclicForwardScatterStatisticalDiffusivePatch, dictionary);


void dsmcCyclicForwardScatterStatisticalDiffusivePatch::readProperties()
{
    Tr0_ = (readScalar(propsDict_.lookup("Tr")));
    velocity_ = (propsDict_.get<vector>("velocity"));
    omegaLim_ = (acos(1 - Foam::sqrt(Foam::sqrt(Tr0_))));
    temperature_ = (readScalar(propsDict_.lookup("temperature")));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcCyclicForwardScatterStatisticalDiffusivePatch::dsmcCyclicForwardScatterStatisticalDiffusivePatch
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
    omegaLim_(acos( 1 - Foam::sqrt(Foam::sqrt(Tr0_)))),
    velocity_(vector(propsDict_.lookup("velocity"))),
    temperature_(readScalar(propsDict_.lookup("temperature"))),
    interp_(x_vals, y_vals, z_vals),
    nReflections_(0),
    nRejections_(0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcCyclicForwardScatterStatisticalDiffusivePatch::~dsmcCyclicForwardScatterStatisticalDiffusivePatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalar dsmcCyclicForwardScatterStatisticalDiffusivePatch::calculateTransmission
(
    const scalar& omega
)
{
    scalar Tr = interp_.interpolate(Tr0_, omega);
    return Tr;
}


void dsmcCyclicForwardScatterStatisticalDiffusivePatch::calculateProperties()
{
    label nReflections = nReflections_;
    label nRejections = nRejections_;

    if(Pstream::parRun())
    {
        reduce(nReflections, sumOp<label>());
        reduce(nRejections, sumOp<label>());
    }

    if(nRejections > 0)
    {
        Info<< "Number of Reflections: " << nReflections
        << ", Number of Rejections: " << nRejections
        << " ratio relfections/(reflections+rejections):"
        << scalar(nReflections)/scalar(nReflections+nRejections)
        << endl;
    }
}

void dsmcCyclicForwardScatterStatisticalDiffusivePatch::initialConfiguration()
{}



void dsmcCyclicForwardScatterStatisticalDiffusivePatch::controlMol
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
//        << " omega = " << omega
//        << " transmission prob.= " << p_
//        << endl;

    Random& rndGen = cloud_.rndGen();

    scalar pRandom = rndGen.sample01<double>();

    if (pRandom >= p_) // reflect
    {
        vector& U = mol.U();

        scalar& ERot = mol.ERot();

        labelList& vibLevel = mol.vibLevel();

        label& ELevel = mol.ELevel();

        label typeId = mol.typeId();

        vector nw = mol.normal();
        nw /= mag(nw);

        // Normal velocity magnitude
        scalar U_dot_nw = U & nw;

        // Wall tangential velocity (flow direction)
        vector Ut = U - U_dot_nw*nw;

        Random& rndGen(cloud_.rndGen());

        while (mag(Ut) < SMALL)
        {
            // If the incident velocity is parallel to the face normal, no
            // tangential direction can be chosen.  Add a perturbation to the
            // incoming velocity and recalculate.

            U = vector
                    (
                            U.x()*(0.8 + 0.2*rndGen.sample01<double>()),
                            U.y()*(0.8 + 0.2*rndGen.sample01<double>()),
                            U.z()*(0.8 + 0.2*rndGen.sample01<double>())
                    );

            U_dot_nw = U & nw;
            Ut = U - U_dot_nw*nw;
        }

        // Wall tangential unit vector
        vector tw1 = Ut/mag(Ut);

        // Other tangential unit vector
        vector tw2 = nw^tw1;

        const scalar& T = temperature_;

        scalar mass = cloud_.constProps(typeId).mass();

        scalar rotationalDof = cloud_.constProps(typeId).rotationalDegreesOfFreedom();

        scalar vibrationalDof = cloud_.constProps(typeId).vibrationalDegreesOfFreedom();


        U =
                sqrt(physicoChemical::k.value()*T/mass)
                *(
                        rndGen.GaussNormal<double>()*tw1
                        + rndGen.GaussNormal<double>()*tw2
                        - sqrt(-2.0*log(max(1 - rndGen.sample01<double>(), VSMALL)))*nw
                );


        ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);


        vibLevel = cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);


        ELevel = cloud_.equipartitionElectronicLevel
                (
                        T,
                        cloud_.constProps(typeId).degeneracyList(),
                        cloud_.constProps(typeId).electronicEnergyList(),
                        typeId
                );

        U += velocity_;

        mol.U() = U;
        mol.ERot() = ERot;
        mol.vibLevel() = vibLevel;
        mol.ELevel() = ELevel;

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



void dsmcCyclicForwardScatterStatisticalDiffusivePatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


void dsmcCyclicForwardScatterStatisticalDiffusivePatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    readProperties();
}


void dsmcCyclicForwardScatterStatisticalDiffusivePatch::autoBoundaryModelMap
(
    const mapPolyMesh& mapper
)
{

}

bool dsmcCyclicForwardScatterStatisticalDiffusivePatch::reflected()
{
    return reflected_;
}


} // End namespace Foam

// ************************************************************************* //
