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

Measures overall temperature, including vibrational temperature, for a single 
species gas or a gas mixture and writes the results to a volume scalar field 
that can be viewed in Paraview.

Translational, rotatational and vibrational temperature field will also be 
written automatically.

Boundary fields are measured in conjunction with the boundaryMeasurements
class and are also written.

\*---------------------------------------------------------------------------*/

#include "dsmcVolFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcVolFields, 1);

addToRunTimeSelectionTable(dsmcField, dsmcVolFields, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcVolFields::dsmcVolFields
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    n_(),
    t1_(),
    t2_(),
    sampleInterval_(1),
    sampleCounter_(0),
    mfpReferenceTemperature_(273.0),
    fieldName_(propsDict_.lookup("fieldName")),
//    rhoNInstantaneous_
//    (
//        IOobject
//        (
//                "rhoNInstantaneous_"+ fieldName_,
//                time_.time().timeName(),
//                mesh_,
//                IOobject::NO_READ,
//                IOobject::AUTO_WRITE
//        ),
//        mesh_,
//        dimensionedScalar(dimless, 0.0)
//    ),
    dsmcRhoN_
    (
        IOobject
        (
            "dsmcRhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
    dsmcRhoNMean_
    (
        IOobject
        (
            "dsmcRhoNMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, VSMALL)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, VSMALL)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, VSMALL)
    ),
    p_
    (
        IOobject
        (
            "p_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("p", dimPressure, VSMALL)
    ),
    translationalT_
    (
        IOobject
        (
            "translationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, VSMALL)
    ),
    rotationalT_
    (
        IOobject
        (
            "rotationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, VSMALL)
    ),
    vibrationalT_
    (
        IOobject
        (
            "vibrationalT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, VSMALL)
    ),
    electronicT_
    (
        IOobject
        (
            "electronicT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, VSMALL)
    ),
    overallT_
    (
        IOobject
        (
            "overallT_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, VSMALL)
    ),
    q_
    (
        IOobject
        (
            "surfaceHeatTransfer_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(1, 0, -3, 0, 0), VSMALL)
    ),
    tau_
    (
        IOobject
        (
            "surfaceShearStress_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, VSMALL)
    ),
    MFP_
    (
        IOobject
        (
            "variableHardSphereMeanFreePath_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, VSMALL)
    ),
    //MFPdX_
    MCR_
    (
        IOobject
        (
            "meanCollisionRate_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, 0, -1, 0, 0), Zero)
    ),
    MCT_
    (
        IOobject
        (
            "meanCollisionTime_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, 0, 1, 0, 0), Zero)
    ),
    MCTdt_
    (
        IOobject
        (
            "mctTimeStepRatio_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    MCS_
    (
        IOobject
        (
            "meanCollisionSeparation"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, Zero)
    ),
    SOF_
    (
        IOobject
        (
            "separationOfFreePaths_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, Zero)
    ),
    Ma_
    (
        IOobject
        (
            "Ma_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    cIDist_
    (
        IOobject
        (
            "classIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    cIIDist_
    (
        IOobject
        (
            "classIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    cIIIDist_
    (
        IOobject
        (
            "classIIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    densityError_
    (
        IOobject
        (
            "densityError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    velocityError_
    (
        IOobject
        (
            "velocityError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    temperatureError_
    (
        IOobject
        (
            "temperatureError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    pressureError_
    (
        IOobject
        (
            "pressureError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    dXMFP_
    (
        IOobject
        (
            "cellMeanFreePathRatio_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, vector::zero)
    ),
    UMean_
    (
        IOobject
        (
            "UMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimLength/dimTime, vector::zero)
    ),
    fD_
    (
        IOobject
        (
            "fD_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVector_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            dimensionSet(1, 0, -3, 0, 0),
            vector::zero
        )
    ),
    pressureTensor_
    (
        IOobject
        (
            "pressureTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            dimPressure,
            tensor::zero
        )
    ),
    shearStressTensor_
    (
        IOobject
        (
            "shearStressTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            dimPressure,
            tensor::zero
        )
    ),
    nTimeSteps_(0.0),
    typeIds_(),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNInstantaneous_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoNMeanInt_(mesh_.nCells(), 0.0),
    molsElec_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    rotationalEMean_(mesh_.nCells(), 0.0),
    rotationalDofMean_(mesh_.nCells(), 0.0),
    muu_(mesh_.nCells(), 0.0),
    muv_(mesh_.nCells(), 0.0),
    muw_(mesh_.nCells(), 0.0),
    mvv_(mesh_.nCells(), 0.0),
    mvw_(mesh_.nCells(), 0.0),
    mww_(mesh_.nCells(), 0.0),
    mcc_(mesh_.nCells(), 0.0),
    mccu_(mesh_.nCells(), 0.0),
    mccv_(mesh_.nCells(), 0.0),
    mccw_(mesh_.nCells(), 0.0),
    eu_(mesh_.nCells(), 0.0),
    ev_(mesh_.nCells(), 0.0),
    ew_(mesh_.nCells(), 0.0),
    e_(mesh_.nCells(), 0.0),
    totalvDof_(mesh_.nCells(), 0.0),
    nClassI_(mesh_.nCells(), 0.0),
    nClassII_(mesh_.nCells(), 0.0),
    nClassIII_(mesh_.nCells(), 0.0),
    collisionSeparation_(mesh_.nCells(), 0.0),
    nColls_(mesh_.nCells(), 0.0),
    momentumMean_(mesh.nCells(), vector::zero),
    momentumMeanXnParticle_(mesh.nCells(), vector::zero),
    boundaryCells_(),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    vibT_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
    vDof_(),
    mfp_(),
    mcr_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    qBF_(),
    vibTxvDofBF_(),
    totalvDofBF_(),
    speciesRhoNIntBF_(),
    speciesRhoNElecBF_(),
    momentumBF_(),
    fDBF_(),
    vibrationalEBF_(),
    electronicEBF_(),
    speciesRhoNBF_(),
    mccSpeciesBF_(),
    vibTBF_(),
    vDofBF_(),
    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false),
    densityOnly_(false),
    measureHeatFluxShearStress_(false)
{

    // standard to reading typeIds ------------ 
    const List<word>& molecules (propsDict_.lookup("typeIds"));

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(moleculesReduced.find(moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(cloud_.typeIdList().find(moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcVolFields::dsmcVolFields()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
       
    // ---------------------------------------------------
    
    // Note; outer list is typeIds, inner list is number of cells on the 
    // mesh
    
    vibT_.setSize(typeIds_.size());
        
    forAll(vibT_, i)
    {
        vibT_[i].setSize(mesh_.nCells());
    }
    
    nGroundElectronicLevel_.setSize(typeIds_.size());
        
    forAll(nGroundElectronicLevel_, i)
    {
        nGroundElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
    }
    
    nFirstElectronicLevel_.setSize(typeIds_.size());
        
    forAll(nFirstElectronicLevel_, i)
    {
        nFirstElectronicLevel_[i].setSize(mesh_.nCells(), 0.0);
    }

    vDof_.setSize(typeIds_.size());
    
    forAll(vDof_, i)
    {
        vDof_[i].setSize(mesh_.nCells());
    }
    
    vibrationalETotal_.setSize(typeIds_.size());
    
//     forAll(vibrationalETotal_, i)
//     {
//         vibrationalETotal_[i].setSize(mesh_.nCells());
//     }
//     
    electronicETotal_.setSize(typeIds_.size());
    
    forAll(electronicETotal_, i)
    {
        electronicETotal_[i].setSize(mesh_.nCells(), 0.0);
    }
    
    nParcels_.setSize(typeIds_.size());
    
    forAll(nParcels_, i)
    {
        nParcels_[i].setSize(mesh_.nCells());
    }
    
    nParcelsXnParticle_.setSize(typeIds_.size());
    
    forAll(nParcelsXnParticle_, i)
    {
        nParcelsXnParticle_[i].setSize(mesh_.nCells());
    }
    
    mccSpecies_.setSize(typeIds_.size());
    
    forAll(nParcels_, i)
    {
        mccSpecies_[i].setSize(mesh_.nCells());
    }
    
    mfp_.setSize(typeIds_.size());
    
    forAll(mfp_, i)
    {
        mfp_[i].setSize(mesh_.nCells());
    }
    
    mcr_.setSize(typeIds_.size());
    
    forAll(mcr_, i)
    {
        mcr_[i].setSize(mesh_.nCells());
    }
    
    boundaryCells_.setSize(mesh.boundaryMesh().size());
    
//     const polyPatch& patch = mesh.boundaryMesh()[p];
    
    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh.boundaryMesh()[p];
        
        boundaryCells_[p].setSize(patch.size());
        
        forAll(boundaryCells_[p], c)
        {
            boundaryCells_[p][c] = patch.faceCells()[c];
        }
    }
        
    // initialisation
    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
    linearKEBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());
    rotationalEBF_.setSize(mesh_.boundaryMesh().size());
    rotationalDofBF_.setSize(mesh_.boundaryMesh().size());
    qBF_.setSize(mesh_.boundaryMesh().size());
    fDBF_.setSize(mesh_.boundaryMesh().size());
    vibTxvDofBF_.setSize(mesh_.boundaryMesh().size());
    totalvDofBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNIntBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNElecBF_.setSize(mesh_.boundaryMesh().size());
    
    n_.setSize(mesh_.boundaryMesh().size());
    t1_.setSize(mesh_.boundaryMesh().size());
    t2_.setSize(mesh_.boundaryMesh().size());
        
    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        
        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), vector::zero);
        rotationalEBF_[j].setSize(patch.size(), 0.0);
        rotationalDofBF_[j].setSize(patch.size(), 0.0);
        qBF_[j].setSize(patch.size(), 0.0);
        fDBF_[j].setSize(patch.size(), vector::zero);
        vibTxvDofBF_[j].setSize(patch.size(), 0.0);
        totalvDofBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNIntBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNElecBF_[j].setSize(patch.size(), 0.0);
        
        n_[j].setSize(patch.size(), vector::zero);
        t1_[j].setSize(patch.size(), vector::zero);
        t2_[j].setSize(patch.size(), vector::zero);
    }
    
    calculateWallUnitVectors();
    
    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    speciesRhoNBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());
    vibTBF_.setSize(typeIds_.size());
    vDofBF_.setSize(typeIds_.size());
    
    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());
        
        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            
            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
            vibTBF_[i][j].setSize(patch.size(), 0.0);
            vDofBF_[i][j].setSize(patch.size(), 0.0);
        }
    }
    
    if (propsDict_.found("sampleInterval"))
    {
        sampleInterval_ = 
                        readLabel(propsDict_.lookup("sampleInterval"));
    }    
    if (propsDict_.found("measureClassifications"))
    {
        measureClassifications_ = 
            Switch(propsDict_.lookup("measureClassifications"));
    }
    
    if (propsDict_.found("measureErrors"))
    {
        measureErrors_ = Switch(propsDict_.lookup("measureErrors"));
    }
    
    if (propsDict_.found("densityOnly"))
    {
        densityOnly_ = Switch(propsDict_.lookup("densityOnly"));
    }
    
    if (propsDict_.found("measureHeatFluxShearStress"))
    {
        measureHeatFluxShearStress_ = 
            Switch(propsDict_.lookup("measureHeatFluxShearStress"));
    }
    
    if(propsDict_.found("measureMeanFreePath"))
    {
        measureMeanFreePath_ = Switch(propsDict_.lookup("measureMeanFreePath"));
    }
    
    if(measureMeanFreePath_)
    {
        mfpReferenceTemperature_ = 
            readScalar(propsDict_.lookup("mfpReferenceTemperature"));
//        measureCellVolume_ = true;
    }
    
    
    if (propsDict_.found("averagingAcrossManyRuns"))
    {
        averagingAcrossManyRuns_ = 
            Switch(propsDict_.lookup("averagingAcrossManyRuns"));
        
        // read in stored data from dictionary
        if(averagingAcrossManyRuns_)
        {
            Info << nl << "Averaging across many runs initiated." << nl << endl;

            readIn();
        }         
    }

//    if (propsDict_.found("measureCellVolume"))
//    {
//        measureCellVolume_ =
//            Switch(propsDict_.lookup("measureCellVolume"));
//    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVolFields::~dsmcVolFields()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcVolFields::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "volFieldsMethod_"+fieldName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("rhoNMean", rhoNMean_);
    dict.readIfPresent("rhoMMean", rhoMMean_);
    dict.readIfPresent("linearKEMean", linearKEMean_);
    dict.readIfPresent("momentumMean", momentumMean_);
    dict.readIfPresent("rotationalEMean", rotationalEMean_);
    dict.readIfPresent("rotationalDofMean", rotationalDofMean_);
    dict.readIfPresent("vibrationalETotal", vibrationalETotal_);
    dict.readIfPresent("nParcels", nParcels_);     
    dict.readIfPresent("rhoNMeanInt", rhoNMeanInt_);     
    
    dict.readIfPresent("nTimeSteps", nTimeSteps_);
}

void dsmcVolFields::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "volFieldsMethod_"+fieldName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("rhoNMean", rhoNMean_);
        dict.add("rhoMMean", rhoMMean_);
        dict.add("linearKEMean", linearKEMean_);
        dict.add("momentumMean", momentumMean_);
        dict.add("rotationalEMean", rotationalEMean_);
        dict.add("rotationalDofMean", rotationalDofMean_);
        dict.add("vibrationalETotal", vibrationalETotal_);
        dict.add("nParcels", nParcels_);     
        dict.add("rhoNMeanInt", rhoNMeanInt_);     
        
        dict.add("nTimeSteps", nTimeSteps_); 
        
//        IOstream::streamFormat fmt = time_.time().writeFormat();
//        IOstream::versionNumber ver = time_.time().writeVersion();
//        IOstream::compressionType cmp = time_.time().writeCompression();

        IOstreamOption::streamFormat fmt = time_.time().writeFormat();
        IOstreamOption::versionNumber ver = time_.time().writeVersion();
        IOstreamOption::compressionType cmp = time_.time().writeCompression();

        //dict.regIOobject::writeObject(fmt, ver, cmp);
        dict.regIOobject::writeObject(fmt, ver, cmp, true);

    }
}

void dsmcVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];
       
        if(isA<wallPolyPatch>(pPatch))
        {           
            const vectorField& fC = pPatch.faceCentres();
           
            forAll(n_[patchi], facei)
            {
                n_[patchi][facei] = pPatch.faceAreas()[facei]
                    /mag(pPatch.faceAreas()[facei]);
               
                //- Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                t1_[patchi][facei] = fC[facei] - 
                            mesh_.points()[mesh_.faces()[pPatch.start() 
                            + facei][0]];
                t1_[patchi][facei] /= mag(t1_[patchi][facei]);
               
                //- Other tangential unit vector.  Rescaling in case face is not
                //  flat and n and t1 aren't perfectly orthogonal
                t2_[patchi][facei] = n_[patchi][facei]^t1_[patchi][facei];
                t2_[patchi][facei] /= mag(t2_[patchi][facei]);
            }
        }
    }
}

//- initial conditions
void dsmcVolFields::createField()
{
    Info << "Initialising dsmcVolFields field" << endl;
    
    forAll(vibrationalETotal_, i)
    {
        vibrationalETotal_[i].setSize(cloud_.constProps(typeIds_[i]).
        vibrationalDegreesOfFreedom());
        
        forAll(vibrationalETotal_[i], j)
        {
            vibrationalETotal_[i][j].setSize(mesh_.nCells(),0.0);
        }
    }
}


void dsmcVolFields::autoMap
(
        const mapPolyMesh &mapper
)
{
    const labelList &cellMap = mapper.cellMap();
    const labelList &reverseCellMap = mapper.reverseCellMap();

    InfoInFunction
            << "Auto Map dsmcVolFields." << endl;
    Pout << "   resize fields from "  << mapper.nOldCells() << " to " << mesh_.nCells() <<  endl;

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
    // init
    // resize
    scalarField correctedRhoNMeanInt_;
    scalarField correctedRhoNMean_;
    scalarField correctedRhoNInstantaneous_;
    scalarField correctedRhoMMean_;

    scalarField correctedLinearKEMean_;
    vectorField correctedMomentumMean_;
    scalarField correctedRotationalEMean_;
    scalarField correctedRotationalDofMean_;

    scalarField correctedRhoNMeanXnParticle_;
    scalarField correctedRhoMMeanXnParticle_;
    vectorField correctedMomentumMeanXnParticle_;
    scalarField correctedLinearKEMeanXnParticle_;

    scalarField correctedMuu_;
    scalarField correctedMuv_;
    scalarField correctedMuw_;
    scalarField correctedMvv_;
    scalarField correctedMvw_;
    scalarField correctedMww_;

    scalarField correctedMcc_;
    scalarField correctedMccu_;
    scalarField correctedMccv_;
    scalarField correctedMccw_;

    scalarField correctedTotalvDof_;

    scalarField correctedEu_;
    scalarField correctedEv_;
    scalarField correctedEw_;
    scalarField correctedE_;

    scalarField correctedNClassI_;
    scalarField correctedNClassII_;
    scalarField correctedNClassIII_;

    scalarField correctedCollisionSeparation_;
    scalarField correctedNColls_;
    scalarField correctedMolsElec_;

    List <scalarField> correctedElectronicETotal_; // extensive
    List <scalarField> correctedNParcels_; // extensive
    List <scalarField> correctedMccSpecies_; // extensive
    List <scalarField> correctedNParcelsXnParticle_; // extensive
    List <scalarField> correctedNGroundElectronicLevel_; // intensive approx.
    List <scalarField> correctedNFirstElectronicLevel_; // intensive approx.
    List <scalarField> correctedMfp_; // intensive
    List <scalarField> correctedMcr_; // intensive

    // resize
    correctedRhoNMeanInt_.setSize(mesh_.nCells());
    correctedRhoNMean_.setSize(mesh_.nCells());
    correctedRhoNInstantaneous_.setSize(mesh_.nCells());
    correctedRhoMMean_.setSize(mesh_.nCells());

    correctedLinearKEMean_.setSize(mesh_.nCells());
    correctedMomentumMean_.setSize(mesh_.nCells());
    correctedRotationalEMean_.setSize(mesh_.nCells());
    correctedRotationalDofMean_.setSize(mesh_.nCells());

    correctedRhoNMeanXnParticle_.setSize(mesh_.nCells());
    correctedRhoMMeanXnParticle_.setSize(mesh_.nCells());
    correctedMomentumMeanXnParticle_.setSize(mesh_.nCells());
    correctedLinearKEMeanXnParticle_.setSize(mesh_.nCells());

    correctedMuu_.setSize(mesh_.nCells());
    correctedMuv_.setSize(mesh_.nCells());
    correctedMuw_.setSize(mesh_.nCells());
    correctedMvv_.setSize(mesh_.nCells());
    correctedMvw_.setSize(mesh_.nCells());
    correctedMww_.setSize(mesh_.nCells());

    correctedMcc_.setSize(mesh_.nCells());
    correctedMccu_.setSize(mesh_.nCells());
    correctedMccv_.setSize(mesh_.nCells());
    correctedMccw_.setSize(mesh_.nCells());

    correctedTotalvDof_.setSize(mesh_.nCells());

    correctedEu_.setSize(mesh_.nCells());
    correctedEv_.setSize(mesh_.nCells());
    correctedEw_.setSize(mesh_.nCells());
    correctedE_.setSize(mesh_.nCells());

    correctedNClassI_.setSize(mesh_.nCells());
    correctedNClassII_.setSize(mesh_.nCells());
    correctedNClassIII_.setSize(mesh_.nCells());

    correctedCollisionSeparation_.setSize(mesh_.nCells());
    correctedNColls_.setSize(mesh_.nCells());
    correctedMolsElec_.setSize(mesh_.nCells());

    createField(); // TODO mapping

    correctedElectronicETotal_.setSize(typeIds_.size()); // extensive
    correctedNParcels_.setSize(typeIds_.size()); // extensive
    correctedMccSpecies_.setSize(typeIds_.size()); // extensive
    correctedNParcelsXnParticle_.setSize(typeIds_.size()); // extensive
    correctedNGroundElectronicLevel_.setSize(typeIds_.size()); // intensive approx.
    correctedNFirstElectronicLevel_.setSize(typeIds_.size()); // intensive approx.
    correctedMfp_.setSize(typeIds_.size()); // intensive
    correctedMcr_.setSize(typeIds_.size()); // intensive

    forAll(typeIds_, iD)
    {
        correctedElectronicETotal_[iD].setSize(mesh_.nCells());
        correctedNParcels_[iD].setSize(mesh_.nCells());
        correctedMccSpecies_[iD].setSize(mesh_.nCells());
        correctedNParcelsXnParticle_[iD].setSize(mesh_.nCells());
        correctedNGroundElectronicLevel_[iD].setSize(mesh_.nCells());
        correctedNFirstElectronicLevel_[iD].setSize(mesh_.nCells());
        correctedMfp_[iD].setSize(mesh_.nCells());
        correctedMcr_[iD].setSize(mesh_.nCells());
    }

    const auto &cellsFromCells = mapper.cellsFromCellsMap();
    Pout << "cellsFromCells = " << cellsFromCells << endl;

    // Correct any split cells
    forAll(cellMap, celli)
    {
        const label oldCelli = cellMap[celli];
        const scalar fRatio = mesh_.V()[celli] / mesh_.V0()[oldCelli];
        Pout << "Cell ratio V/V0 = " << fRatio << endl;

//        if (oldCelli >= 0 && nSubCells[oldCelli] == 8)
        if (oldCelli >= 0)
        {
            Pout << "Handle refined cell: celli=" << celli << " oldCelli=" << oldCelli << " nSubCells=" << nSubCells[oldCelli] << endl;
            // Found split cell. Use current volume instead of mapped
            // old one
            correctedRhoNMeanInt_[celli] = rhoNMeanInt_[oldCelli] * fRatio;
            correctedRhoNMean_[celli] = rhoNMean_[oldCelli] * fRatio; // extensive; is numbers per cell
            correctedRhoNInstantaneous_[celli] = rhoNInstantaneous_[oldCelli] * fRatio;
            correctedRhoMMean_[celli] = rhoMMean_[oldCelli] * fRatio; // extensive; is mass per cell

            correctedLinearKEMean_[celli] = linearKEMean_[oldCelli] * fRatio;
            correctedMomentumMean_[celli] = momentumMean_[oldCelli] * fRatio;
            correctedRotationalEMean_[celli] = rotationalEMean_[oldCelli] * fRatio;
            correctedRotationalDofMean_[celli] = rotationalDofMean_[oldCelli] * fRatio;

            correctedRhoNMeanXnParticle_[celli] = rhoNMeanXnParticle_[oldCelli] * fRatio;
            correctedRhoMMeanXnParticle_[celli] = rhoMMeanXnParticle_[oldCelli] * fRatio;
            correctedMomentumMeanXnParticle_[celli] = momentumMeanXnParticle_[oldCelli] * fRatio;
            correctedLinearKEMeanXnParticle_[celli] = linearKEMeanXnParticle_[oldCelli] * fRatio;

            correctedMuu_[celli] = muu_[oldCelli] * fRatio;
            correctedMuv_[celli] = muv_[oldCelli] * fRatio;
            correctedMuw_[celli] = muw_[oldCelli] * fRatio;
            correctedMvv_[celli] = mvv_[oldCelli] * fRatio;
            correctedMvw_[celli] = mvw_[oldCelli] * fRatio;
            correctedMww_[celli] = mww_[oldCelli] * fRatio;
            correctedMcc_[celli] = mcc_[oldCelli] * fRatio;
            correctedMccu_[celli] = mccu_[oldCelli] * fRatio;
            correctedMccv_[celli] = mccv_[oldCelli] * fRatio;
            correctedMccw_[celli] = mccw_[oldCelli] * fRatio;

            correctedTotalvDof_[celli] = totalvDof_[oldCelli] * fRatio;

            correctedEu_[celli] = eu_[oldCelli] * fRatio;
            correctedEv_[celli] = ev_[oldCelli] * fRatio;
            correctedEw_[celli] = ew_[oldCelli] * fRatio;
            correctedE_[celli] = e_[oldCelli] * fRatio;

            correctedNClassI_[celli] = nClassI_[oldCelli] * fRatio;
            correctedNClassII_[celli] = nClassII_[oldCelli] * fRatio;
            correctedNClassIII_[celli] = nClassIII_[oldCelli] * fRatio;

            correctedCollisionSeparation_[celli] = collisionSeparation_[oldCelli]; // intensive
            correctedNColls_[celli] = nColls_[oldCelli] * fRatio;
            correctedMolsElec_[celli] = molsElec_[oldCelli] * fRatio;

            forAll(typeIds_, iD)
            {
                correctedElectronicETotal_[iD][celli] = electronicETotal_[iD][oldCelli] * fRatio;
                correctedNParcels_[iD][celli] = nParcels_[iD][oldCelli] * fRatio;
                correctedMccSpecies_[iD][celli] = mccSpecies_[iD][oldCelli] * fRatio;
                correctedNParcelsXnParticle_[iD][celli] = nParcelsXnParticle_[iD][oldCelli] * fRatio;
                correctedNGroundElectronicLevel_[iD][celli] = nGroundElectronicLevel_[iD][oldCelli]; // intensive
                correctedNFirstElectronicLevel_[iD][celli] = nFirstElectronicLevel_[iD][oldCelli]; // intensive
                correctedMfp_[iD][celli] = mfp_[iD][oldCelli]; // intensive
                correctedMcr_[iD][celli] = mcr_[iD][oldCelli]; // intensive
            }
        }
    }

    for (const auto &s: cellsFromCells)
    {
        // Reset unsplit cell
        const label celli = s.index();
        // sum up extensive vars
        const labelList &oldCells = s.masterObjects();
        Pout << "Handle unrefined cells. celli=" << celli << " oldCells=" << oldCells << endl;
        for (const label oldCelli: oldCells)
        {

            correctedRhoNMeanInt_[celli] += rhoNMeanInt_[oldCelli];
            correctedRhoNMean_[celli] = rhoNMean_[oldCelli]; // intensive
            correctedRhoNInstantaneous_[celli] += rhoNInstantaneous_[oldCelli];
            correctedRhoMMean_[celli] = rhoMMean_[oldCelli]; // intensive

            correctedLinearKEMean_[celli] += linearKEMean_[oldCelli];
            correctedMomentumMean_[celli] += momentumMean_[oldCelli];
            correctedRotationalEMean_[celli] += rotationalEMean_[oldCelli];
            correctedRotationalDofMean_[celli] += rotationalDofMean_[oldCelli];

            correctedRhoNMeanXnParticle_[celli] += rhoNMeanXnParticle_[oldCelli];
            correctedRhoMMeanXnParticle_[celli] += rhoMMeanXnParticle_[oldCelli];
            correctedMomentumMeanXnParticle_[celli] += momentumMeanXnParticle_[oldCelli];
            correctedLinearKEMeanXnParticle_[celli] += linearKEMeanXnParticle_[oldCelli];

            correctedMuu_[celli] += muu_[oldCelli];
            correctedMuv_[celli] += muv_[oldCelli];
            correctedMuw_[celli] += muw_[oldCelli];
            correctedMvv_[celli] += mvv_[oldCelli];
            correctedMvw_[celli] += mvw_[oldCelli];
            correctedMww_[celli] += mww_[oldCelli];
            correctedMcc_[celli] += mcc_[oldCelli];
            correctedMccu_[celli] += mccu_[oldCelli];
            correctedMccv_[celli] += mccv_[oldCelli];
            correctedMccw_[celli] += mccw_[oldCelli];

            correctedTotalvDof_[celli] += totalvDof_[oldCelli];

            correctedEu_[celli] += eu_[oldCelli];
            correctedEv_[celli] += ev_[oldCelli];
            correctedEw_[celli] += ew_[oldCelli];
            correctedE_[celli] += e_[oldCelli];

            correctedNClassI_[celli] += nClassI_[oldCelli];
            correctedNClassII_[celli] += nClassII_[oldCelli];
            correctedNClassIII_[celli] += nClassIII_[oldCelli];

            correctedCollisionSeparation_[celli] = collisionSeparation_[oldCelli]; // intensive
            correctedNColls_[celli] += nColls_[oldCelli];
            correctedMolsElec_[celli] += molsElec_[oldCelli];

            forAll(typeIds_, iD)
            {
                correctedElectronicETotal_[iD][celli] += electronicETotal_[iD][oldCelli];
                correctedNParcels_[iD][celli] += nParcels_[iD][oldCelli];
                correctedMccSpecies_[iD][celli] += mccSpecies_[iD][oldCelli];
                correctedNParcelsXnParticle_[iD][celli] += nParcelsXnParticle_[iD][oldCelli];
                correctedNGroundElectronicLevel_[iD][celli] += nGroundElectronicLevel_[iD][oldCelli]; // intensive
                correctedNFirstElectronicLevel_[iD][celli] += nFirstElectronicLevel_[iD][oldCelli]; // intensive
                correctedMfp_[iD][celli] = mfp_[iD][oldCelli]; // intensive
                correctedMcr_[iD][celli] = mcr_[iD][oldCelli]; // intensive
            }

        }
    }

    rhoNMeanInt_ = correctedRhoNMeanInt_;
    rhoNMean_ = correctedRhoNMean_;
//    rhoNInstantaneous_.field() = correctedRhoNInstantaneous_;
    rhoNInstantaneous_ = correctedRhoNInstantaneous_;
    rhoMMean_ = correctedRhoMMean_;

    linearKEMean_ = correctedLinearKEMean_;
    momentumMean_ = correctedMomentumMean_;
    rotationalEMean_ = correctedRotationalEMean_;
    rotationalDofMean_ = correctedRotationalDofMean_;

    rhoNMeanXnParticle_ = correctedRhoNMeanXnParticle_;
    rhoMMeanXnParticle_ = correctedRhoMMeanXnParticle_;
    momentumMeanXnParticle_ = correctedMomentumMeanXnParticle_;
    linearKEMeanXnParticle_ = correctedLinearKEMeanXnParticle_;

    muu_ = correctedMuu_;
    muv_ = correctedMuv_;
    muw_ = correctedMuw_;
    mvv_ = correctedMvv_;
    mvw_ = correctedMvw_;
    mww_ = correctedMww_;
    mcc_ = correctedMcc_;
    mccu_ = correctedMccu_;
    mccv_ = correctedMccv_;
    mccw_ = correctedMccw_;

    totalvDof_ = correctedTotalvDof_;

    eu_ = correctedEu_;
    ev_ = correctedEv_;
    ew_ = correctedEw_;
    e_ = correctedE_;

    nClassI_ = correctedNClassI_;
    nClassII_ = correctedNClassII_;
    nClassIII_ = correctedNClassIII_;

    collisionSeparation_ = correctedCollisionSeparation_; // intensive
    nColls_ = correctedNColls_;
    molsElec_ = correctedMolsElec_;

    forAll(typeIds_, iD)
    {
        electronicETotal_[iD] = correctedElectronicETotal_[iD];
        nParcels_[iD] = correctedNParcels_[iD];
        mccSpecies_[iD] = correctedMccSpecies_[iD];
        nParcelsXnParticle_[iD] = correctedNParcelsXnParticle_[iD];
        nGroundElectronicLevel_[iD] = correctedNGroundElectronicLevel_[iD]; // intensive
        nFirstElectronicLevel_[iD] = correctedNFirstElectronicLevel_[iD]; // intensive
        mfp_[iD] = correctedMfp_[iD]; // intensive
        mcr_[iD] = correctedMcr_[iD]; // intensive
    }


    // boundaries
    Pout << "   autoMap boundary elements" << endl;
    Pout << "   resize boundaryMesh" << endl;
    boundaryCells_.setSize(mesh_.boundaryMesh().size());

    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
    linearKEBF_.setSize(mesh_.boundaryMesh().size());
    rotationalEBF_.setSize(mesh_.boundaryMesh().size());
    rotationalDofBF_.setSize(mesh_.boundaryMesh().size());
    qBF_.setSize(mesh_.boundaryMesh().size());
    vibTxvDofBF_.setSize(mesh_.boundaryMesh().size());
    totalvDofBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNIntBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNElecBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());
    fDBF_.setSize(mesh_.boundaryMesh().size());

    n_.setSize(mesh_.boundaryMesh().size());
    t1_.setSize(mesh_.boundaryMesh().size());
    t2_.setSize(mesh_.boundaryMesh().size());

    forAll(rhoNBF_, j)
    {
        const polyPatch &patch = mesh_.boundaryMesh()[j];
        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        rotationalEBF_[j].setSize(patch.size(), 0.0);
        rotationalDofBF_[j].setSize(patch.size(), 0.0);
        qBF_[j].setSize(patch.size(), 0.0);
        vibTxvDofBF_[j].setSize(patch.size(), 0.0);
        totalvDofBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNIntBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNElecBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), vector::zero);
        fDBF_[j].setSize(patch.size(), vector::zero);

        n_[j].setSize(patch.size(), vector::zero);
        t1_[j].setSize(patch.size(), vector::zero);
        t2_[j].setSize(patch.size(), vector::zero);
    }

    calculateWallUnitVectors();

    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    speciesRhoNBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());
    vibTBF_.setSize(typeIds_.size());
    vDofBF_.setSize(typeIds_.size());

    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());

        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch &patch = mesh_.boundaryMesh()[j];

            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
            vibTBF_[i][j].setSize(patch.size(), 0.0);
            vDofBF_[i][j].setSize(patch.size(), 0.0);
        }
    }
    //
}
    /*
    // obtain boundary measurements
    Pout << "   reobtain boundary measurements" << endl;

    forAll(cloud_.boundaryFluxMeasurements().rhoNBF(), i)
    {
        const label& iD = typeIds_.find(i);

        forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i], j)
        {
            forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j], k)
            {
                if(iD != -1)
                {
                    rhoNBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                    rhoMBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k];
                    linearKEBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k];
                    momentumBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k];
                    rotationalEBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k];
                    rotationalDofBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k];
                    qBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().qBF()[i][j][k];
                    fDBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().fDBF()[i][j][k];
                    speciesRhoNBF_[iD][j][k] +=
                            cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                    vibrationalEBF_[iD][j][k] +=
                            cloud_.boundaryFluxMeasurements().vibrationalEBF()[i][j][k];
                    electronicEBF_[iD][j][k] +=
                            cloud_.boundaryFluxMeasurements().electronicEBF()[i][j][k];
                    mccSpeciesBF_[iD][j][k] +=
                            cloud_.boundaryFluxMeasurements().mccSpeciesBF()[i][j][k];
                    speciesRhoNIntBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().rhoNIntBF()[i][j][k];
                    speciesRhoNElecBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements().rhoNElecBF()[i][j][k];
                }
            }
        }
    }


}
*/
//
//
//    // boundaries
//    Pout << "   autoMap boundary elements" << endl;
//    Pout << "   resize boundaryMesh" << endl;
//    boundaryCells_.setSize(mesh_.boundaryMesh().size());
//
//    rhoNBF_.setSize(mesh_.boundaryMesh().size());
//    rhoMBF_.setSize(mesh_.boundaryMesh().size());
//    linearKEBF_.setSize(mesh_.boundaryMesh().size());
//    rotationalEBF_.setSize(mesh_.boundaryMesh().size());
//    rotationalDofBF_.setSize(mesh_.boundaryMesh().size());
//    qBF_.setSize(mesh_.boundaryMesh().size());
//    vibTxvDofBF_.setSize(mesh_.boundaryMesh().size());
//    totalvDofBF_.setSize(mesh_.boundaryMesh().size());
//    speciesRhoNIntBF_.setSize(mesh_.boundaryMesh().size());
//    speciesRhoNElecBF_.setSize(mesh_.boundaryMesh().size());
//    momentumBF_.setSize(mesh_.boundaryMesh().size());
//    fDBF_.setSize(mesh_.boundaryMesh().size());
//
//    n_.setSize(mesh_.boundaryMesh().size());
//    t1_.setSize(mesh_.boundaryMesh().size());
//    t2_.setSize(mesh_.boundaryMesh().size());
//
//    forAll(rhoNBF_, j)
//    {
//        const polyPatch& patch = mesh_.boundaryMesh()[j];
//        rhoNBF_[j].setSize(patch.size(),0.0);
//        rhoMBF_[j].setSize(patch.size(),0.0);
//        linearKEBF_[j].setSize(patch.size(),0.0);
//        rotationalEBF_[j].setSize(patch.size(),0.0);
//        rotationalDofBF_[j].setSize(patch.size(),0.0);
//        qBF_[j].setSize(patch.size(),0.0);
//        vibTxvDofBF_[j].setSize(patch.size(),0.0);
//        totalvDofBF_[j].setSize(patch.size(),0.0);
//        speciesRhoNIntBF_[j].setSize(patch.size(),0.0);
//        speciesRhoNElecBF_[j].setSize(patch.size(),0.0);
//        momentumBF_[j].setSize(patch.size(),vector::zero);
//        fDBF_[j].setSize(patch.size(),vector::zero);
//
//        n_[j].setSize(patch.size(), vector::zero);
//        t1_[j].setSize(patch.size(), vector::zero);
//        t2_[j].setSize(patch.size(), vector::zero);
//    }
//
//    calculateWallUnitVectors();
//
//    vibrationalEBF_.setSize(typeIds_.size());
//    electronicEBF_.setSize(typeIds_.size());
//    speciesRhoNBF_.setSize(typeIds_.size());
//    mccSpeciesBF_.setSize(typeIds_.size());
//    vibTBF_.setSize(typeIds_.size());
//    vDofBF_.setSize(typeIds_.size());
//
//    forAll(vibrationalEBF_, i)
//    {
//        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
//        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
//        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
//        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
//        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
//        vDofBF_[i].setSize(mesh_.boundaryMesh().size());
//
//        forAll(vibrationalEBF_[i], j)
//        {
//            const polyPatch& patch = mesh_.boundaryMesh()[j];
//
//            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
//            electronicEBF_[i][j].setSize(patch.size(), 0.0);
//            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
//            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
//            vibTBF_[i][j].setSize(patch.size(), 0.0);
//            vDofBF_[i][j].setSize(patch.size(), 0.0);
//        }
//    }
//    /*
//    // obtain boundary measurements
//    Pout << "   reobtain boundary measurements" << endl;
//
//    forAll(cloud_.boundaryFluxMeasurements().rhoNBF(), i)
//    {
//        const label& iD = typeIds_.find(i);
//
//        forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i], j)
//        {
//            forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j], k)
//            {
//                if(iD != -1)
//                {
//                    rhoNBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
//                    rhoMBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k];
//                    linearKEBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k];
//                    momentumBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k];
//                    rotationalEBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k];
//                    rotationalDofBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k];
//                    qBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().qBF()[i][j][k];
//                    fDBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().fDBF()[i][j][k];
//                    speciesRhoNBF_[iD][j][k] +=
//                            cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
//                    vibrationalEBF_[iD][j][k] +=
//                            cloud_.boundaryFluxMeasurements().vibrationalEBF()[i][j][k];
//                    electronicEBF_[iD][j][k] +=
//                            cloud_.boundaryFluxMeasurements().electronicEBF()[i][j][k];
//                    mccSpeciesBF_[iD][j][k] +=
//                            cloud_.boundaryFluxMeasurements().mccSpeciesBF()[i][j][k];
//                    speciesRhoNIntBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().rhoNIntBF()[i][j][k];
//                    speciesRhoNElecBF_[j][k] +=
//                            cloud_.boundaryFluxMeasurements().rhoNElecBF()[i][j][k];
//                }
//            }
//        }
//    }




void dsmcVolFields::calculateField()
{
//    Info << "dsmcVolFields::calculateField()" << nl << endl;
    sampleCounter_++;
    
    const scalar& nParticle = cloud_.nParticle();

    rhoNInstantaneous_ = 0.0;
//    Info << "rhoNInstantaneous_.field() = " << rhoNInstantaneous_.field() << endl;
    //rhoNInstantaneous_.setSize(mesh_.nCells(), 0);
    
    if(sampleInterval_ <= sampleCounter_)
    {
        nTimeSteps_ += 1.0;

        
        if(densityOnly_)
        {
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                label iD = typeIds_.find(p.typeId());

                if(iD != -1)
                {
                    const label& cell = p.cell();
                    const scalar& mass = cloud_.constProps(p.typeId()).mass();

                    rhoNMean_[cell] += 1.0;
                    rhoNInstantaneous_[cell] += 1.0;
                    
                    if(cloud_.axisymmetric())
                    {
                        const point& cC = cloud_.mesh().cellCentres()[cell];
                        
                        scalar radius = cC.y();
                        
                        scalar RWF = 1.0;
                        
                        RWF = 1.0 + 
                            cloud_.maxRWF()*(radius/cloud_.radialExtent());
                        
                        rhoNMeanXnParticle_[cell] += (RWF*nParticle);
                        rhoMMeanXnParticle_[cell] += (mass*RWF*nParticle);
                    }
                    else
                    {
                        rhoNMeanXnParticle_[cell] += nParticle;
                        rhoMMeanXnParticle_[cell] += (mass*nParticle);
                    }

                }
            }
        }
        else
        {
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                const label& iD = typeIds_.find(p.typeId());
//                Pout << "iD: " << iD << endl;

                if(iD != -1)
                {
                    const label& cell = p.cell();
                    const scalar& mass = cloud_.constProps(p.typeId()).mass();
                    const scalar& massByMagUsq = mass*sqr(mag(p.U()));
                    const scalarList& electronicEnergies = 
                        cloud_.constProps(typeIds_[iD]).electronicEnergyList();
                    const scalar& rotationalDof = 
                    cloud_.constProps(p.typeId()).rotationalDegreesOfFreedom();
                    const scalar& xVel = p.U().x();
                    const scalar& yVel = p.U().y();
                    const scalar& zVel = p.U().z();
                    const vector& U = p.U();


                    scalarList EVib
                    (
                        cloud_.constProps(typeIds_[iD])
                        .vibrationalDegreesOfFreedom()
                    );
                    
                    if(EVib.size() > 0)
                    {
                        forAll(EVib, i)
                        {
                        
                            EVib[i] = p.vibLevel()[i]
                                     *physicoChemical::k.value()                
                                     *cloud_.constProps(p.typeId()).thetaV()[i];
                            
                            vibrationalETotal_[iD][i][cell] +=           
                                        p.vibLevel()[i]
                                        *physicoChemical::k.value()
                                     *cloud_.constProps(p.typeId()).thetaV()[i];
                        }
                    }
                    rhoNMean_[cell] += 1.0;
                    rhoNInstantaneous_[cell] += 1.0;
                    rhoMMean_[cell] += mass;
                    linearKEMean_[cell] += mass*(U & U);
                    momentumMean_[cell] += mass*U;
                    rotationalEMean_[cell] += p.ERot();
                    rotationalDofMean_[cell] += rotationalDof;
                    electronicETotal_[iD][cell] += 
                                        electronicEnergies[p.ELevel()];
                    nParcels_[iD][cell] += 1.0;
                    mccSpecies_[iD][cell] += massByMagUsq;

                    
                    if(cloud_.axisymmetric())
                    {
                        const point& cC = cloud_.mesh().cellCentres()[cell];
                        
                        scalar radius = cC.y();
                        
                        scalar RWF = 1.0;
                        
                        RWF = 1.0 + 
                            cloud_.maxRWF()*(radius/cloud_.radialExtent());
                        
                        nParcelsXnParticle_[iD][cell] += (RWF*nParticle);
                        rhoNMeanXnParticle_[cell] += (RWF*nParticle);
                        rhoMMeanXnParticle_[cell] += (mass*RWF*nParticle);
                        momentumMeanXnParticle_[cell] += 
                                    (mass*(U)*RWF*nParticle);
                        linearKEMeanXnParticle_[cell] += (mass*(U & 
                                                        U)*RWF*nParticle);
                    }
                    else
                    {

                        nParcelsXnParticle_[iD][cell] += nParticle;
                        rhoNMeanXnParticle_[cell] += nParticle;
                        rhoMMeanXnParticle_[cell] += (mass*nParticle);
                        momentumMeanXnParticle_[cell] += 
                                                (mass*(U)*nParticle);
                        linearKEMeanXnParticle_[cell] += (mass*(U & 
                                                        U)*nParticle);
                    }

//                    Info << "specific" << endl;
                    muu_[cell] += mass*sqr(xVel);
                    muv_[cell] += mass*( xVel * yVel );
                    muw_[cell] += mass*( xVel * zVel );
                    mvv_[cell] += mass*sqr(yVel);
                    mvw_[cell] += mass*( yVel * zVel );
                    mww_[cell] += mass*sqr(zVel);
                    
                    mcc_[cell] += massByMagUsq;
                    mccu_[cell] += massByMagUsq*(xVel);
                    mccv_[cell] += massByMagUsq*(yVel);
                    mccw_[cell] += massByMagUsq*(xVel);

                    
                    scalar vibEn = 0.0;
                    
                    if(EVib.size() > 0)
                    {
                        forAll(EVib, v)
                        {
                            vibEn += EVib[v];
                        }
                    }


                    eu_[cell] += ( p.ERot() + vibEn )*(xVel);
                    ev_[cell] += ( p.ERot() + vibEn )*(yVel);
                    ew_[cell] += ( p.ERot() + vibEn )*(zVel);
                    e_[cell] += ( p.ERot() + vibEn );

                    if(rotationalDof > VSMALL)
                    {
                        rhoNMeanInt_[cell] += 1.0;
                    }
                    
                    label nElecLevels =                 
                       cloud_.constProps(p.typeId()).numberOfElectronicLevels();
                    
                    if(nElecLevels > 1)
                    {
                        molsElec_[cell] += 1.0;
                        
                        if(p.ELevel() == 0)
                        {
                            nGroundElectronicLevel_[iD][cell]++;
                        }
                        if(p.ELevel() == 1)
                        {
                            nFirstElectronicLevel_[iD][cell]++;
                        }
                    }
                    if(measureClassifications_)
                    {

                        label classification = p.classification();
                        
                        if(classification == 0)
                        {
                            nClassI_[cell] += 1.0;
                        }
                        
                        if(classification == 1)
                        {
                            nClassII_[cell] += 1.0;
                        }
                        
                        if(classification == 2)
                        {
                            nClassIII_[cell] += 1.0;
                        }
                    }
                }
            }
//            Pout << "linearKEMeanXnParticle_ = " << linearKEMeanXnParticle_ << endl;
//            Pout << "momentumMeanXnParticle_ = " << momentumMeanXnParticle_ << endl;
//
//            Pout << "rhoMMeanXnParticle_ = " <<rhoMMeanXnParticle_ << endl;
//            Pout << "rhoNMeanXnParticle_ = " <<rhoNMeanXnParticle_ << nl << endl;
//
//            Pout << "rhoNMean_ = " << rhoNMean_ << endl;
//            Pout << "rhoNInstantaneous_ = " << rhoNInstantaneous_ << nl << endl;
//            ^ disabled: dumps the whole per-cell field on every rank every time
//              step, which makes parallel runs I/O bound (GB-sized logs)

            
            // obtain collision quality measurements
            Info << "Obtain collision quality measurements" << endl;
            forAll(cloud_.cellPropMeasurements().collisionSeparation(), cell)
            {
                collisionSeparation_[cell] += 
                    cloud_.cellPropMeasurements().collisionSeparation()[cell];
                nColls_[cell] += cloud_.cellPropMeasurements().nColls()[cell];
            }
            
            // obtain boundary measurements
            Info << "Obtain boundary measurements" << endl;
            
            forAll(cloud_.boundaryFluxMeasurements().rhoNBF(), i)
            {
                const label& iD = typeIds_.find(i);
            
                forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i], j)
                {                
                    forAll(cloud_.boundaryFluxMeasurements().rhoNBF()[i][j], k)
                    {
                        if(iD != -1)
                        { 
                            rhoNBF_[j][k] += 
                            cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                            rhoMBF_[j][k] += 
                            cloud_.boundaryFluxMeasurements().rhoMBF()[i][j][k];
                            linearKEBF_[j][k] += 
                        cloud_.boundaryFluxMeasurements().linearKEBF()[i][j][k];
                            momentumBF_[j][k] += 
                        cloud_.boundaryFluxMeasurements().momentumBF()[i][j][k];
                            rotationalEBF_[j][k] += 
                    cloud_.boundaryFluxMeasurements().rotationalEBF()[i][j][k];
                            rotationalDofBF_[j][k] += 
                cloud_.boundaryFluxMeasurements().rotationalDofBF()[i][j][k];
                            qBF_[j][k] += 
                               cloud_.boundaryFluxMeasurements().qBF()[i][j][k];
                            fDBF_[j][k] += 
                              cloud_.boundaryFluxMeasurements().fDBF()[i][j][k];
                            speciesRhoNBF_[iD][j][k] += 
                            cloud_.boundaryFluxMeasurements().rhoNBF()[i][j][k];
                            vibrationalEBF_[iD][j][k] += 
                    cloud_.boundaryFluxMeasurements().vibrationalEBF()[i][j][k];
                            electronicEBF_[iD][j][k] += 
                    cloud_.boundaryFluxMeasurements().electronicEBF()[i][j][k];
                            mccSpeciesBF_[iD][j][k] += 
                    cloud_.boundaryFluxMeasurements().mccSpeciesBF()[i][j][k];
                            speciesRhoNIntBF_[j][k] += 
                        cloud_.boundaryFluxMeasurements().rhoNIntBF()[i][j][k];
                            speciesRhoNElecBF_[j][k] += 
                        cloud_.boundaryFluxMeasurements().rhoNElecBF()[i][j][k];
                        }
                    }
                }
            }
        }
        sampleCounter_ = 0;

    }
    
    if(time_.time().outputTime())
    {
        Info << "Mark fields for writing to time directory." << nl << endl;
        const scalar& nAvTimeSteps = nTimeSteps_;
        
        if(densityOnly_)
        {
            forAll(rhoNMean_, cell)
            {
                if(rhoNMean_[cell] > VSMALL)
                {
                    const scalar& cellVolume = mesh_.cellVolumes()[cell];
                    
                    dsmcRhoNMean_[cell] = rhoNMean_[cell]/(nAvTimeSteps);
                    
                    rhoN_[cell] = 
                        rhoNMeanXnParticle_[cell]/(nAvTimeSteps*cellVolume);
                    
                    rhoM_[cell] = 
                        rhoMMeanXnParticle_[cell]/(nAvTimeSteps*cellVolume);

                }
                else
                {
                     // not zero so that weighted decomposition still works
                     WarningInFunction
                        << "rhoNMean_[cell] < VSMALL: Set variable "
                           "(dsmcRhoNMean_ rhoN_ rhoM_) to non-zero value 0.001" << nl << endl;
                    dsmcRhoNMean_[cell] = 0.001;
                    rhoN_[cell] = 0.001;
                    rhoM_[cell] = 0.001;
                }
                
                if(rhoNInstantaneous_[cell] > VSMALL)
                {
                    dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
                    dsmcRhoN_[cell] = 0.001;
                }
            }
        }
        else
        {                  
            scalarField vibT(mesh_.nCells(), scalar(0.0));
            scalarField vibTForOverallT(mesh_.nCells(), scalar(0.0));
            scalarField Cp(mesh_.nCells(), scalar(0.0));
            scalarField Cv(mesh_.nCells(), scalar(0.0));
            scalarField molecularMass(mesh_.nCells(), scalar(0.0));
            scalarField Cv_p(mesh_.nCells(), scalar(0.0));
            scalarField totalvDof(mesh_.nCells(), scalar(0.0));
            scalarField totalvDofOverall(mesh_.nCells(), scalar(0.0));

            UMean_.writeOpt(IOobject::AUTO_WRITE);
            translationalT_.writeOpt(IOobject::AUTO_WRITE);
            rotationalT_.writeOpt(IOobject::AUTO_WRITE);
            vibrationalT_.writeOpt(IOobject::AUTO_WRITE);
            p_.writeOpt(IOobject::AUTO_WRITE);
            overallT_.writeOpt(IOobject::AUTO_WRITE);
            electronicT_.writeOpt(IOobject::AUTO_WRITE);
            Ma_.writeOpt(IOobject::AUTO_WRITE);
            fD_.writeOpt(IOobject::AUTO_WRITE);
            q_.writeOpt(IOobject::AUTO_WRITE);
            tau_.writeOpt(IOobject::AUTO_WRITE);

//            rhoNInstantaneous_.writeOpt(IOobject::AUTO_WRITE);
//            Pout <<  "measureHeatFluxShearStress_" << endl;
            if(measureHeatFluxShearStress_)
            {
                pressureTensor_.writeOpt(IOobject::AUTO_WRITE);
                shearStressTensor_.writeOpt(IOobject::AUTO_WRITE);
                heatFluxVector_.writeOpt(IOobject::AUTO_WRITE);

            }
//            Pout <<  "measureClassifications_" << endl;
            if(measureClassifications_)
            {
                cIDist_.writeOpt(IOobject::AUTO_WRITE);
                cIIDist_.writeOpt(IOobject::AUTO_WRITE);
                cIIIDist_.writeOpt(IOobject::AUTO_WRITE);

            }
//            Pout <<  "measureMeanFreePath_" << endl;
            if(measureMeanFreePath_)
            {
                MFP_.writeOpt(IOobject::AUTO_WRITE);
                dXMFP_.writeOpt(IOobject::AUTO_WRITE);
                MCR_.writeOpt(IOobject::AUTO_WRITE);
                MCT_.writeOpt(IOobject::AUTO_WRITE);
                MCTdt_.writeOpt(IOobject::AUTO_WRITE);
                MCS_.writeOpt(IOobject::AUTO_WRITE);
                SOF_.writeOpt(IOobject::AUTO_WRITE);
            }

//            Pout <<  "measureErrors_" << endl;
            if(measureErrors_)
            {
                densityError_.writeOpt(IOobject::AUTO_WRITE);
                velocityError_.writeOpt(IOobject::AUTO_WRITE);
                temperatureError_.writeOpt(IOobject::AUTO_WRITE);
                pressureError_.writeOpt(IOobject::AUTO_WRITE);
            }



            Pout
                << "    writing mean components" << endl;
            forAll(rhoNMean_, cell)
            {                
                if(rhoNMean_[cell] > VSMALL)
                {                  
                    const scalar& cellVolume = mesh_.cellVolumes()[cell];
//                    Pout << "cellVolume[cell]= " << cellVolume << endl;
//                    Pout << "rhoMMeanXnParticle_[cell]= " << rhoMMeanXnParticle_[cell] << endl;
//                    Pout << "nAvTimeSteps" << nAvTimeSteps << endl;

                    if (rhoNMean_[cell]/(nAvTimeSteps) <= SMALL) {
                        Pout << "dsmcRhoNMean_[cell]= " << rhoNMean_[cell] / (nAvTimeSteps) << endl;
                    }
                    dsmcRhoNMean_[cell] = rhoNMean_[cell]/(nAvTimeSteps);
                    if ((rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume) <= SMALL) {
                        Pout << "rhoN_[cell]= " << (rhoNMeanXnParticle_[cell]) / (nAvTimeSteps * cellVolume) << endl;
                    }
                    rhoN_[cell] = (rhoNMeanXnParticle_[cell]) / (nAvTimeSteps * cellVolume);

                    if ((rhoMMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume) <= SMALL) {
                        Pout << "rhoM_[cell]= " << (rhoMMeanXnParticle_[cell]) / (nAvTimeSteps * cellVolume) << endl;
                    }
                    rhoM_[cell] =
                            rhoMMeanXnParticle_[cell] / (nAvTimeSteps * cellVolume);

                    if (rhoMMeanXnParticle_[cell] / (cellVolume * nAvTimeSteps) <= SMALL) {
                        Pout << "rhoMMean_[cell]= " << rhoMMeanXnParticle_[cell] / (cellVolume * nAvTimeSteps) << endl;
                    }
                    scalar rhoMMean =
                                rhoMMeanXnParticle_[cell] / (cellVolume * nAvTimeSteps);

                    if (Foam::mag(momentumMeanXnParticle_[cell]) /
                        (rhoMMean * cellVolume * nAvTimeSteps) <= SMALL) {
                        Pout << "UMean_[cell]= " << momentumMeanXnParticle_[cell] /
                                                    (rhoMMean * cellVolume * nAvTimeSteps) << endl;
                    }
                    UMean_[cell] = momentumMeanXnParticle_[cell] /
                                   (rhoMMean * cellVolume * nAvTimeSteps);

                    if (0.5 * linearKEMeanXnParticle_[cell]
                        / (cellVolume * nAvTimeSteps)<= SMALL) {
                        Pout << "linearKEMean_[cell]= " << 0.5 * linearKEMeanXnParticle_[cell]
                                                           / (cellVolume * nAvTimeSteps) << endl;
                    }
                    scalar linearKEMean = 0.5 * linearKEMeanXnParticle_[cell]
                                              / (cellVolume * nAvTimeSteps);


                    if (rhoNMeanXnParticle_[cell] / (cellVolume * nAvTimeSteps) <= SMALL) {
                        Pout << "rhoNMean[cell]= " << rhoNMeanXnParticle_[cell] / (cellVolume * nAvTimeSteps) << endl;
                    }
                    scalar rhoNMean =
                            rhoNMeanXnParticle_[cell] / (cellVolume * nAvTimeSteps);


//                    if (2.0 / (3.0 * physicoChemical::k.value() * rhoNMean)
//                        * (linearKEMean - 0.5 * rhoMMean *
//                                          (UMean_[cell] & UMean_[cell])
//                        )<= SMALL) {

                        //linearKEMean = 0.5 * linearKEMeanXnParticle_[cell]
                        //               / (cellVolume * nAvTimeSteps)
//                    Pout << "linearKEMeanXnParticle_[cell=" << cell << "]= "<<  linearKEMeanXnParticle_[cell] << endl;
//                    Pout << "linearKEMean= " << linearKEMean << endl;
//                    Pout << "rhoMMean=" << rhoMMean << rhoMMeanXnParticle_[cell]/ (cellVolume * nAvTimeSteps) << endl;
//                    Pout << "rhoMMeanXnParticle_[cell=" << cell << "] = " << rhoMMeanXnParticle_[cell] << endl;
//                    Pout << "momentumMeanXnParticle_[cell=" << cell << "] = " << momentumMeanXnParticle_[cell] << endl;
//                    Pout << "UMean[cell=" << cell << "] = " <<  UMean_[cell] << endl;
//                    Pout << "(0.5 * rhoMMean * (UMean_[cell] & UMean_[cell]) = "
//                            << (0.5 * rhoMMean * (UMean_[cell] & UMean_[cell])) << endl;
//                    Pout << "translationalT_[cell=" << cell << "]= " << 2.0 / (3.0 * physicoChemical::k.value() * rhoNMean)
//                                                         * (linearKEMean - 0.5 * rhoMMean *
//                                                                           (UMean_[cell] & UMean_[cell])
//                                                         ) << endl;

                    translationalT_[cell] =
                            2.0 / (3.0 * physicoChemical::k.value() * rhoNMean)
                            * (linearKEMean - 0.5 * rhoMMean *
                                              (UMean_[cell] & UMean_[cell])
                            );
                    // Statistical noise in a cell with barely any parcels can
                    // make the fluctuation energy (and hence T and p) slightly
                    // negative. Clamp it: a negative T makes later sqrt() calls
                    // raise FE_INVALID and trap with FOAM_SIGFPE.
                    if (translationalT_[cell] < SMALL)
                    {
                        translationalT_[cell] = SMALL;
                    }
//                    Pout << "p_[cell=" << cell << "]= " << rhoN_[cell] * physicoChemical::k.value()
//                                                           * translationalT_[cell] << nl << endl;

                    if (rhoN_[cell] * physicoChemical::k.value()
                        * translationalT_[cell]<= SMALL) {
                        Pout << "ERROR: p_[cell=" << cell << "]= " << rhoN_[cell] * physicoChemical::k.value()
                                                * translationalT_[cell] << endl;
                    }
                    p_[cell] =
                            rhoN_[cell] * physicoChemical::k.value()
                            * translationalT_[cell];
                }
                else
                {
                    WarningInFunction
                        << "Small value (" << rhoNMean_[cell]
                        << ") of rhoNMean_ in cell " << cell
                        << "." << endl;
                    // not zero so that weighted decomposition still works
                    rhoNMean_[cell] = 0.001;
                    dsmcRhoNMean_[cell] = 0.001;
                    rhoN_[cell] = VSMALL;
                    rhoM_[cell] = VSMALL;
                    UMean_[cell] = vector::zero;
                    translationalT_[cell] = 0.0;           
                    p_[cell] = 0.0;
                }


                if(rhoNInstantaneous_[cell] > VSMALL)
                {
                    dsmcRhoN_[cell] = rhoNInstantaneous_[cell];
                }
                else
                {
//                    Pout << "rhoNInstantaneous_[cell=" << cell << "] = " << rhoNInstantaneous_[cell] << endl;
                    dsmcRhoN_[cell] = 0.001;
                }
//                Pout << "rotationalDofMean_" << endl;
                if(rotationalDofMean_[cell] > VSMALL && nAvTimeSteps > VSMALL)
                {
                    scalar rotationalEMean = rotationalEMean_[cell] / 
                                                nAvTimeSteps;
                    scalar rotationalDofMean = rotationalDofMean_[cell] / 
                                                nAvTimeSteps;

                    rotationalT_[cell] =  (2.0/physicoChemical::k.value())
                                           *(rotationalEMean/rotationalDofMean);
                }
                else
                {
                    rotationalT_[cell] = 0.0;
                }
               
                scalarList degreesOfFreedomSpecies(typeIds_.size(),0.0);
                scalarList vibTID(vibrationalETotal_.size(), 0.0);
                
                List<scalarList> dofMode;
                List<scalarList> vibTMode;
                
                dofMode.setSize(typeIds_.size());
                vibTMode.setSize(typeIds_.size());
                
                forAll(dofMode, iD)
                { 
                    dofMode[iD].setSize(cloud_.constProps(typeIds_[iD]).
                                        vibrationalDegreesOfFreedom(), 0.0);
                    
                    vibTMode[iD].setSize(cloud_.constProps(typeIds_[iD]).
                                        vibrationalDegreesOfFreedom(), 0.0);
                }
//                Pout << "vibrationalETotal_" << endl;
                forAll(vibrationalETotal_, iD)
                {
                    forAll(vibrationalETotal_[iD], v)
                    {
                        if(vibrationalETotal_[iD][v][cell] > VSMALL
                                    && nParcels_[iD][cell] > VSMALL
                                    && dofMode.size() > VSMALL)
                        {        
                            scalar thetaV = 
                                cloud_.constProps(typeIds_[iD]).thetaV()[v];
                            
                            scalar vibrationalEMean = 
                                    vibrationalETotal_[iD][v][cell]
                                    /nParcels_[iD][cell];
                            
                            scalar iMean = 
                                    vibrationalEMean
                                    /(physicoChemical::k.value()*thetaV);
                            
                            vibTMode[iD][v] = thetaV / log(1.0 + (1.0/iMean));

                            dofMode[iD][v] = 
                                (2.0*thetaV/vibTMode[iD][v]) 
                                / (exp(thetaV/vibTMode[iD][v]) - 1.0);
                        }
                    }
                    
                    forAll(dofMode[iD], v)
                    {
                        degreesOfFreedomSpecies[iD] += 
                                    dofMode[iD][v];
                    }
                    
                    forAll(dofMode[iD], v)
                    {
                        if(degreesOfFreedomSpecies[iD] > VSMALL)
                        {
                            vibTID[iD] += 
                                vibTMode[iD][v]
                                *(dofMode[iD][v]
                                /degreesOfFreedomSpecies[iD]);
                        }
                    }

                    
                    totalvDof[cell] += degreesOfFreedomSpecies[iD];
                    
                    if(rhoNMeanInt_[cell] > VSMALL 
                        && rhoNMean_[cell] > VSMALL 
                        && nParcels_[iD][cell] > VSMALL)
                    {
                        scalar fraction = 
                                nParcels_[iD][cell]
                                /rhoNMeanInt_[cell];
                        
                        scalar fractionOverall = 
                                nParcels_[iD][cell]
                                /rhoNMean_[cell];
                        
                        totalvDofOverall[cell] += 
                                totalvDof[cell]
                                *(fractionOverall/fraction);
                        
                        vibT[cell] += vibTID[iD]*fraction;
                    }
                }

                vibrationalT_[cell] = vibT[cell];
                
//                 forAll(vibrationalETotal_, iD)
//                 {
//                     if(vibrationalETotal_[iD][cell] > VSMALL && 
//nParcels_[iD][cell] > VSMALL)
//                     {        
//                         const scalar& thetaV = 
//cloud_.constProps(typeIds_[iD]).thetaV();
//                         
//                         scalar vibrationalEMean = 
//(vibrationalETotal_[iD][cell]/nParcels_[iD][cell]);
//                         
//                         scalar iMean = 
//vibrationalEMean/(physicoChemical::k.value()*thetaV);
//                         
//                         vibT_[iD][cell] = thetaV / log(1.0 + (1.0/iMean));
//                         
//                         scalar fraction = 
//nParcels_[iD][cell]/rhoNMeanInt_[cell];
//                         
//                         vibT[cell] += vibT_[iD][cell]*fraction;
//                         
//                         vDof_[iD][cell] = 
//fraction*(2.0*thetaV/vibT_[iD][cell]) / (exp(thetaV/vibT_[iD][cell]) - 1.0);
//                         
//                         totalvDof_[cell] += vDof_[iD][cell];
//                     }
//                 }
// 
//                 vibrationalT_[cell] = vibT[cell];
                
                // electronic temperature
                scalar totalEDof = 0.0;
                scalar elecT = 0.0;

//                Pout << "nParcels_" << endl;
                forAll(nParcels_, iD)
                {
                    const scalarList& electronicEnergies = 
                        cloud_.constProps(typeIds_[iD]).electronicEnergyList();
                    const labelList& degeneracies = 
                        cloud_.constProps(typeIds_[iD]).degeneracyList();

                    if(nGroundElectronicLevel_[iD][cell] > VSMALL && 
                        nFirstElectronicLevel_[iD][cell] > VSMALL && 
                        nFirstElectronicLevel_[iD][cell]*degeneracies[0] != 
                        nGroundElectronicLevel_[iD][cell]*degeneracies[1])
                    {
                        
                        scalar elecTID = 
                            (electronicEnergies[1]-electronicEnergies[0])/
                            (
                                physicoChemical::k.value()*
                                log((nGroundElectronicLevel_[iD][cell]*         
                                 degeneracies[1])/
                                (nFirstElectronicLevel_[iD][cell]*
                                degeneracies[0]))   
                            );

                    
                        scalar fraction = nParcels_[iD][cell]/molsElec_[cell];
                            
                        if(elecTID > VSMALL)
                        {
                            elecT += fraction*elecTID;
                        }
                        
                        
                        scalar eDof = (2.0*(electronicETotal_[iD][cell]/
                             nParcels_[iD][cell]))/
                             (physicoChemical::k.value()*elecTID);

                        totalEDof += fraction*eDof;
                    }
                    
    //                 label nElectronicLevels = 
//cloud_.constProps(typeIds_[iD]).numberOfElectronicLevels();
    //                 
    //                 if(nElectronicLevels > 1 && nParcels_[iD][cell] > VSMALL 
//&& molsElec_[cell] > VSMALL)
    //                 {
    //                     const scalarList& electronicEnergies = 
//cloud_.constProps(typeIds_[iD]).electronicEnergyList();
    //                     const labelList& degeneracies = 
//cloud_.constProps(typeIds_[iD]).degeneracyList();
    //                     
    //                     scalar speciesTransT = 
//(1.0/(3.0*physicoChemical::k.value()))
    //                                             *(
    //                                                 
//(mccSpecies_[iD][cell]/(nParcels_[iD][cell]))
    //                                                 - (
    //                                                     
//cloud_.constProps(typeIds_[iD]).mass()*mag(UMean_[cell])*mag(UMean_[cell])
    //                                                 )
    //                                             );
    //                     
    //                     scalar fraction = 
//nParcels_[iD][cell]/molsElec_[cell];
    //                     
    //                     if(speciesTransT > SMALL && 
//electronicETotal_[iD][cell] > VSMALL)
    //                     {
    //                         scalar sum1 = 0.0;
    //                         scalar sum2 = 0.0;
    //                         
    //                         forAll(electronicEnergies, ii)
    //                         {
    //                             sum1 += 
// degeneracies[ii]*exp(-electronicEnergies[ii]/(physicoChemical::k.value()*
// speciesTransT));
    //                             sum2 += 
// degeneracies[ii]*(electronicEnergies[ii]/(physicoChemical::k.value()*
// speciesTransT))
    //                                         
// *exp(-electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT));
    //                         }
    //                         
    //                         if(sum2 > VSMALL && sum1 > VSMALL)
    //                         {
    //                             scalar elecTID = 
// 
// (electronicETotal_[iD][cell]/(physicoChemical::k.value()*nParcels_[iD][cell]))*( 
// sum1/sum2);
    //                             
    //                             if(elecTID > SMALL && elecTID < GREAT)
    //                             {
    //                                 elecT += fraction*elecTID;
    //                                 
    //                                 scalar eDof = 
// (2.0*(electronicETotal_[iD][cell]/nParcels_[iD][cell]))/(physicoChemical::k.
// value()*speciesTransT);
    //                                 
    //                                 totalEDof += fraction*eDof;
    //                             }
    //                         }
    //                     }
    //                 }
                }

                electronicT_[cell] = elecT;

                scalar nRotDof = 0.0;

//                Pout << "rhoNMean_[cell] > VSMALL?" << endl;
                
                if(rhoNMean_[cell] > VSMALL)
                {
                    nRotDof = rotationalDofMean_[cell] / rhoNMean_[cell];
                }
                
                overallT_[cell] = ( 
                                    (3.0*translationalT_[cell]) 
                                    + (nRotDof*rotationalT_[cell]) 
                                    + (totalvDof_[cell]*vibrationalT_[cell])
                                    + (totalEDof*electronicT_[cell])
                                   ) /
                                    (3.0 + nRotDof + totalvDof_[cell] + 
                                        totalEDof);
//                Pout << "measureHeatFluxShearStress_" << endl;
                if(measureHeatFluxShearStress_)
                {                    
                    if(rhoNMean_[cell] > VSMALL)
                    {
                        pressureTensor_[cell].xx() = rhoN_[cell]*
                        ( 
                            muu_[cell]/(rhoNMean_[cell]) - 
                            (
                                (rhoMMean_[cell]/(rhoNMean_[cell]))
                                *UMean_[cell].x()*UMean_[cell].x()
                            )
                        );
                        pressureTensor_[cell].xy() = rhoN_[cell]*
                        ( 
                            muv_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell])))
                            *UMean_[cell].x()*UMean_[cell].y() 
                            
                        );
                        pressureTensor_[cell].xz() = rhoN_[cell]*
                        ( 
                            muw_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell]))
                            *UMean_[cell].x()*UMean_[cell].z())  
                        );
                        pressureTensor_[cell].yx() = pressureTensor_[cell].xy();
                        pressureTensor_[cell].yy() = rhoN_[cell]*
                        ( 
                            mvv_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell])))
                            *UMean_[cell].y()*UMean_[cell].y()
                        );
                        pressureTensor_[cell].yz() = rhoN_[cell]*
                        ( 
                            mvw_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell]))
                            *UMean_[cell].y()*UMean_[cell].z())
                        );
                        pressureTensor_[cell].zx() = pressureTensor_[cell].xz();
                        pressureTensor_[cell].zy() = pressureTensor_[cell].yz();
                        pressureTensor_[cell].zz() = rhoN_[cell]*
                        (
                            mww_[cell]/(rhoNMean_[cell]) - 
                            ((rhoMMean_[cell]/(rhoNMean_[cell]))
                            *UMean_[cell].z()*UMean_[cell].z()) 
                        );
                                                
                        scalar scalarPressure = (1.0/3.0)*
                                                (pressureTensor_[cell].xx() + 
                                                pressureTensor_[cell].yy() + 
                                                pressureTensor_[cell].zz());
                        
                        shearStressTensor_[cell] = -pressureTensor_[cell];
                        shearStressTensor_[cell].xx() += scalarPressure;
                        shearStressTensor_[cell].yy() += scalarPressure;
                        shearStressTensor_[cell].zz() += scalarPressure;
                        
                       //terms involving pressure tensor should not be 
                       //multiplied by the number density (see Bird corrigendum)
                        
                        heatFluxVector_[cell].x() = rhoN_[cell]*
                        (                     
                            0.5*(mccu_[cell]/(rhoNMean_[cell])) - 
                            0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                            UMean_[cell].x() + eu_[cell]/(rhoNMean_[cell]) - 
                            (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].x()
                        ) - 
                            pressureTensor_[cell].xx()*UMean_[cell].x() - 
                            pressureTensor_[cell].xy()*UMean_[cell].y() - 
                            pressureTensor_[cell].xz()*UMean_[cell].z();
                                                
                        heatFluxVector_[cell].y() = rhoN_[cell]*
                        (                          
                            0.5*(mccv_[cell]/(rhoNMean_[cell])) - 
                            0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                            UMean_[cell].y() + ev_[cell]/(rhoNMean_[cell])- 
                            (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].y()
                        ) - 
                            pressureTensor_[cell].yx()*UMean_[cell].x() - 
                            pressureTensor_[cell].yy()*UMean_[cell].y() - 
                            pressureTensor_[cell].yz()*UMean_[cell].z();
                                                
                        heatFluxVector_[cell].z() = rhoN_[cell]*
                        (                  
                            0.5*(mccw_[cell]/(rhoNMean_[cell])) - 
                            0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                            UMean_[cell].z() + ew_[cell]/(rhoNMean_[cell]) - 
                            (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].z()
                        ) - 
                            pressureTensor_[cell].zx()*UMean_[cell].x() - 
                            pressureTensor_[cell].zy()*UMean_[cell].y() - 
                            pressureTensor_[cell].zz()*UMean_[cell].z();
                    }
                    else
                    {
                        pressureTensor_[cell] = tensor::zero;
                        shearStressTensor_[cell] = tensor::zero;
                        heatFluxVector_[cell] = vector::zero;
                    }
                }

                totalvDof_[cell] = 0.0;
                 
                forAll(nParcels_, iD)  
                {
                    label typeId = typeIds_[iD];

                    if(rhoNMean_[cell] > VSMALL)
                    {
                        molecularMass[cell] += 
                                        cloud_.constProps(typeId).mass()*
                                        (nParcels_[iD][cell]/rhoNMean_[cell]);
                        Cp[cell] += (5.0 + 
                         cloud_.constProps(typeId).rotationalDegreesOfFreedom())
                            *(nParcels_[iD][cell]/rhoNMean_[cell]);
                        Cv[cell] += (3.0 + 
                        cloud_.constProps(typeId).rotationalDegreesOfFreedom())
                            *(nParcels_[iD][cell]/rhoNMean_[cell]);
                    }
                }
                
                Cv_p[cell] = Cv[cell]/6.02214129e23;

                scalar gasConstant = 0.0;
                scalar gamma = 0.0;
                scalar speedOfSound = 0.0;
                
                if(molecularMass[cell] > VSMALL)
                {
                    gasConstant = // R = k/m
                        physicoChemical::k.value()/molecularMass[cell]; 
                }
                
                if(Cv[cell] > VSMALL)
                {
                    gamma = Cp[cell]/Cv[cell]; // gamma = cP/cV
                }
                
                if(gamma > VSMALL && gasConstant > VSMALL 
                    && translationalT_[cell] > VSMALL)
                {
                    speedOfSound = 
                        sqrt(gamma*gasConstant*translationalT_[cell]);
                }
                
                if(speedOfSound > VSMALL)
                {
                    Ma_[cell] = mag(UMean_[cell])/speedOfSound;
                }
                else
                {
                    Ma_[cell] = 0.0;
                }
//                Pout << "measureMeanFreePath_::4004" << endl;
                if(measureMeanFreePath_)
                {
                    forAll(mfp_, iD)
                    {
                        label qspec = 0;
                        
                        for (qspec=0; qspec<typeIds_.size(); qspec++)
                        {
                            scalar dPQ = 
                                0.5*(cloud_.constProps(typeIds_[iD]).d() + 
                                cloud_.constProps(typeIds_[qspec]).d());
                            scalar omegaPQ = 
                                0.5*(cloud_.constProps(typeIds_[iD]).omega() + 
                                cloud_.constProps(typeIds_[qspec]).omega());
                            scalar massRatio = 
                                cloud_.constProps(typeIds_[iD]).mass()/
                                cloud_.constProps(typeIds_[qspec]).mass();
                            
                            if(nParcels_[qspec][cell] > VSMALL && 
                                        translationalT_[cell] > VSMALL)
                            {
                                scalar nDensQ = 
                                    (nParcelsXnParticle_[qspec][cell])/
                                    (mesh_.cellVolumes()[cell]*nTimeSteps_);
                                    
                                scalar reducedMass = 
                                (
                                    cloud_.constProps(typeIds_[iD]).mass()*
                                    cloud_.constProps(typeIds_[qspec]).mass()  
                                )/ 
                                (
                                    cloud_.constProps(typeIds_[iD]).mass()+
                                    cloud_.constProps(typeIds_[qspec]).mass() 
                                );
                                
                                mfp_[iD][cell] += 
                                    (
                                        pi*dPQ*dPQ*nDensQ*
                                        pow(mfpReferenceTemperature_/
                                        translationalT_[cell], omegaPQ-0.5)
                                        *sqrt(1.0+massRatio)
                                    ); //Bird, eq (4.76)
                                
                                mcr_[iD][cell] += 
                                (
                                    2.0*sqrt(pi)*dPQ*dPQ*nDensQ*
                                    pow(translationalT_[cell]/
                                    mfpReferenceTemperature_, 1.0-omegaPQ)  
                                    *sqrt(2.0*physicoChemical::k.value()*
                                    mfpReferenceTemperature_/reducedMass)
                                ); // //Bird, eq (4.74)
                            }
                        }
                        
                        if(mfp_[iD][cell] > VSMALL)
                        {
                            mfp_[iD][cell] = 1.0/mfp_[iD][cell];
                        }
                    }

                    MFP_[cell] = 0.0;
                    //MFPdX_[cell] = 0.0;
                    dXMFP_[cell] = vector::zero;
                    MCR_[cell] = 0.0;
                    MCT_[cell] = 0.0;
                    MCTdt_[cell] = 0.0;
                    MCS_[cell] = 0.0;
                    
                    if(nColls_[cell] > VSMALL)
                    {
                        MCS_[cell] = collisionSeparation_[cell]/nColls_[cell];
                    }
                    else
                    {
                       MCS_[cell] = GREAT; 
                    }
                    
                    forAll(mfp_, iD)
                    {
                        if(rhoN_[cell] > VSMALL)
                        {                    
                            scalar nDensP = (nParcelsXnParticle_[iD][cell])/
                                        (mesh_.cellVolumes()[cell]*nTimeSteps_);

                            //Bird, eq (4.77)
                            MFP_[cell] += mfp_[iD][cell]*nDensP/rhoN_[cell];

                            //Bird, eq (1.38)
                            MCR_[cell] += mcr_[iD][cell]*nDensP/rhoN_[cell];
                        }
                    }

/*                    if(MFP_[cell] < VSMALL)
                    {
                        MFP_[cell] = GREAT;
                    }*/
    
                    const scalar& deltaT = mesh_.time().deltaTValue();

                    if(MCR_[cell] > VSMALL)
                    {
                        MCT_[cell] = 1.0/MCR_[cell];
                        MCTdt_[cell] = MCT_[cell]/deltaT;
                    }
                    else
                    {
                        MCT_[cell] = GREAT;
                        MCTdt_[cell] = GREAT;
                    }
                
                    forAll(mfp_, iD)
                    {
                        mfp_[iD][cell] = 0.0;
                        mcr_[iD][cell] = 0.0;
                    }

                    if(MFP_[cell] != GREAT)
                    {
                        scalar largestXCellDimension = 0.0;
                        scalar largestYCellDimension = 0.0;
                        scalar largestZCellDimension = 0.0;

                        const labelList& pLabels(mesh_.cells()[cell].
                                                    labels(mesh_.faces()));
                        pointField pLocal(pLabels.size(), vector::zero);

                        forAll (pLabels, pointi)
                        {
                            pLocal[pointi] = mesh_.points()[pLabels[pointi]];
                        }

                        scalarField dimension;

                        dimension.setSize(3, 0.0);

                        dimension[0] =  Foam::max(pLocal & vector(1,0,0)) -
                                            Foam::min(pLocal & vector(1,0,0));
                        dimension[1] = Foam::max(pLocal & vector(0,1,0)) -
                                            Foam::min(pLocal & vector(0,1,0));
                        dimension[2] = Foam::max(pLocal & vector(0,0,1)) -
                                       Foam::min(pLocal & vector(0,0,1));

                        largestXCellDimension = dimension[0];
                        largestYCellDimension = dimension[1];
                        largestZCellDimension = dimension[2];

                        label dim = 0;

//                        for (dim=0; dim<dimension.size(); dim++)
//                        {
//                            if(dimension[dim] > largestCellDimension)
//                            {
//                                largestCellDimension = dimension[dim];
//                            }
//                        }

                        //MFPdX_[cell] = MFP_[cell]/largestCellDimension;

                        dXMFP_[cell] = vector(dimension[0]/MFP_[cell],
                                              dimension[1]/MFP_[cell],
                                              dimension[2]/MFP_[cell]);


                        if(MFP_[cell] > VSMALL && MCS_[cell] > VSMALL)
                        {
                            SOF_[cell] = MCS_[cell]/MFP_[cell];
                        }
                    }
                    else
                    {
                        //MFPdX_[cell] = GREAT;
                        dXMFP_[cell] = vector(SMALL, SMALL, SMALL);
                        SOF_[cell] = GREAT;
                    }
                }
                
                if(measureClassifications_)
                {            
                    if(rhoNMean_[cell] > VSMALL)
                    {
                        cIDist_[cell] = nClassI_[cell]/rhoNMean_[cell];
                        cIIDist_[cell] = nClassII_[cell]/rhoNMean_[cell];
                        cIIIDist_[cell] = nClassIII_[cell]/rhoNMean_[cell];
                    }
                }
                
                if(measureErrors_)
                {
                    if(dsmcRhoNMean_[cell] > VSMALL && 
                        Ma_[cell] > VSMALL && gamma > VSMALL 
                        && Cv_p[cell] > VSMALL)
                    {
                        densityError_[cell] =1.0/
                                        sqrt(dsmcRhoNMean_[cell]*nTimeSteps_);
                        velocityError_[cell] = 
                            (1.0/sqrt(dsmcRhoNMean_[cell]*nTimeSteps_))
                            *(1.0/(Ma_[cell]*sqrt(gamma)));
                        temperatureError_[cell] = 
                                (1.0/sqrt(dsmcRhoNMean_[cell]*nTimeSteps_))*   
                                    sqrt(physicoChemical::k.value()/Cv_p[cell]);
                        pressureError_[cell] = sqrt(gamma)/
                                        sqrt(dsmcRhoNMean_[cell]*nTimeSteps_);
                    }
                    
                }
            }
            
            List<scalarField> vibTBF(mesh_.boundaryMesh().size());
            List<scalarField> molMassBoundary(mesh_.boundaryMesh().size());
            List<scalarField> CpBoundary(mesh_.boundaryMesh().size());
            List<scalarField> CvBoundary(mesh_.boundaryMesh().size());
            List<scalarField> Cv_pBoundary(mesh_.boundaryMesh().size());

            // computing boundary measurements
            Pout << "computing boundary measurements" << endl;
            forAll(rhoNBF_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                vibTBF[j].setSize(patch.size(), 0.0);
                molMassBoundary[j].setSize(patch.size(), 0.0);
                CpBoundary[j].setSize(patch.size(), 0.0);
                CvBoundary[j].setSize(patch.size(), 0.0); 
                Cv_pBoundary[j].setSize(patch.size(), 0.0);
                
                if(isA<wallPolyPatch>(patch))
                {                               
                    forAll(rhoN_.boundaryField()[j], k)
                    {                        
                        rhoN_.boundaryFieldRef()[j][k] =
                                rhoNBF_[j][k]*nParticle/nAvTimeSteps; // replaced .boundaryFied() with .boundaryFieldRef()
                        rhoM_.boundaryFieldRef()[j][k] =
                                rhoMBF_[j][k]*nParticle/nAvTimeSteps; // ***
                        
                        if(rhoM_.boundaryField()[j][k] > VSMALL)
                        {
                            UMean_.boundaryFieldRef()[j][k] =
                                    momentumBF_[j][k]*nParticle/
                                    (rhoM_.boundaryField()[j][k]*nAvTimeSteps); // ***
                        }
                        else
                        {
                            UMean_.boundaryFieldRef()[j][k] = vector::zero; // ***
                        }
                            
                        scalar rhoMMean = rhoMBF_[j][k]*nParticle/nAvTimeSteps;
                        scalar linearKEMean = linearKEBF_[j][k]*
                                                nParticle/nAvTimeSteps;
                        scalar rhoNMean = rhoNBF_[j][k]*nParticle/nAvTimeSteps;
                        
                        if(rhoNMean > VSMALL)
                        {
                            translationalT_.boundaryFieldRef()[j][k] =
                                2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                            *(linearKEMean - 0.5*rhoMMean*
                                            (
                                                UMean_.boundaryField()[j][k] & 
                                                UMean_.boundaryField()[j][k])
                                            ); // ***
                        }
                        else
                        {
                            translationalT_.boundaryFieldRef()[j][k] = 0.0; // ***
                        }
                        
                        if(rotationalDofBF_[j][k] > VSMALL)
                        {
                            rotationalT_.boundaryFieldRef()[j][k] =
                                        (2.0/physicoChemical::k.value())
                                *(rotationalEBF_[j][k]/rotationalDofBF_[j][k]);
                        }
                        else
                        {
                            rotationalT_.boundaryFieldRef()[j][k] = 0.0; // ***
                        }
                        
//                         forAll(vibrationalEBF_, i)
//                         {
//                             if(rhoNBF_[j][k] > VSMALL)
//                             {                       
//                                 molMassBoundary[j][k] +=  
//cloud_.constProps(typeIds_[i]).mass()
//                                             
//*(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);
//                                             
//                                 
// CpBoundary[j][k] += (5.0 + 
// cloud_.constProps(typeIds_[i]).rotationalDegreesOfFreedom())
//                                             
// *(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);
//                                             
//                                 
// CvBoundary[j][k] += (3.0 + 
// cloud_.constProps(typeIds_[i]).rotationalDegreesOfFreedom())
//                                             
// *(speciesRhoNBF_[i][j][k]/rhoNBF_[j][k]);
//                             }
//                             
//                             if(vibrationalEBF_[i][j][k] > VSMALL && 
// speciesRhoNBF_[i][j][k] > VSMALL && speciesRhoNIntBF_[j][k] > VSMALL)
//                             {        
//                                 const scalar& thetaV = 
// cloud_.constProps(typeIds_[i]).thetaV();
//                                 
//                                 scalar vibrationalEMean = 
// (vibrationalEBF_[i][j][k]/speciesRhoNBF_[i][j][k]);
//                                 
//                                 scalar iMean = 
// vibrationalEMean/(physicoChemical::k.value()*thetaV);
//                                 
//                                 vibTBF_[i][j][k] = thetaV / log(1.0 + 
// (1.0/iMean));
//                                 
//                                 scalar fraction = 
// speciesRhoNBF_[i][j][k]/speciesRhoNIntBF_[j][k];
//                                 
//                                 vDofBF_[i][j][k] = 
// fraction*(2.0*thetaV/vibTBF_[i][j][k]) / (exp(thetaV/vibTBF_[i][j][k]) - 
// 1.0);
//                                 
//                                 vibTBF[j][k] += fraction*vibTBF_[i][j][k];
//                                 
//                                 totalvDofBF_[j][k] += vDofBF_[i][j][k];
//                             }
// 
//                         }
//                         
//                         if(totalvDofBF_[j][k] > VSMALL)
//                         {
//                             vibrationalT_.boundaryField()[j][k] = 
// vibTBF[j][k];
//                         }
//                         else
//                         {
//                             vibrationalT_.boundaryField()[j][k] = 0.0;
//                         }
                        
                        // electronic temperature
                        scalar totalEDof = 0.0;
                        scalar elecT = 0.0;
                            
    //                     forAll(electronicEBF_, i)
    //                     {
    //                         label nElectronicLevels = 
// cloud_.constProps(typeIds_[i]).numberOfElectronicLevels();
    //                         
    //                         if(nElectronicLevels > 1 && 
// speciesRhoNBF_[i][j][k] > VSMALL && speciesRhoNElecBF_[j][k] > VSMALL)
    //                         {
    //                             const scalarList& electronicEnergies = 
// cloud_.constProps(typeIds_[i]).electronicEnergyList();
    //                             const labelList& degeneracies = 
// cloud_.constProps(typeIds_[i]).degeneracyList();
    //                             
    //                             scalar speciesTransT = 
// (1.0/(3.0*physicoChemical::k.value()))
    //                                                     *(
    //                                                         
// ((mccSpeciesBF_[i][j][k]/(speciesRhoNBF_[i][j][k])))
    //                                                         - (
    //                                                             
// (cloud_.constProps(typeIds_[i]).mass()
    //                                                             
// )*mag(UMean_.boundaryField()[j][k])*mag(UMean_.boundaryField()[j][k]))
    //                                                     );
    //                             
    //                             scalar fraction = 
// speciesRhoNBF_[i][j][k]/speciesRhoNElecBF_[j][k];
    //                             
    //                             if(speciesTransT > VSMALL)
    //                             {
    //                                 scalar sum1 = 0.0;
    //                                 scalar sum2 = 0.0;
    //                                 
    //                                 forAll(electronicEnergies, ii)
    //                                 {
    //                                     sum1 += 
// degeneracies[ii]*exp(-electronicEnergies[ii]/(physicoChemical::k.value()*
// speciesTransT));
    //                                     sum2 += 
// degeneracies[ii]*(electronicEnergies[ii]/(physicoChemical::k.value()*
// speciesTransT))
    //                                                 
// *exp(-electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT));
    //                                 }
    //                                 
    //                                 if(sum2 > VSMALL && sum1> VSMALL)
    //                                 {
    //                                     scalar elecTID = 
// 
// 
//(electronicEBF_[i][j][k]/(physicoChemical::k.value()*speciesRhoNBF_[i][j][k]))
//* ( // sum1/sum2);
    //                                     
    //                                     if(elecTID > SMALL && elecTID < 
// GREAT)
    //                                     {
    //                                         elecT += fraction*elecTID;
    //                                         
    //                                         scalar eDof = 
// (2.0*(electronicEBF_[i][j][k]/speciesRhoNBF_[i][j][k]))/(physicoChemical::k.
// value()*speciesTransT);
    //                                         
    //                                         totalEDof += fraction*eDof;
    //                                     }
    //                                 }
    //                             }
    //                         }
    //                     }
                        
                        electronicT_.boundaryFieldRef()[j][k] = elecT; // changed .boundaryField() to .boundaryFieldRef()
                        
                        scalar nRotDof = 0.0;
                        
                        if(rhoNBF_[j][k] > VSMALL)
                        {
                            nRotDof = rotationalDofBF_[j][k] / rhoNBF_[j][k];
                        }
                        
                        overallT_.boundaryFieldRef()[j][k] =
                        (                      
                            (3.0*translationalT_.boundaryField()[j][k]) + 
                            (nRotDof*rotationalT_.boundaryField()[j][k]) + 
                                (totalvDofBF_[j][k]*
                                vibrationalT_.boundaryField()[j][k]) + 
                            (totalEDof*elecT)
                        ) /
                        (3.0 + nRotDof + totalvDofBF_[j][k]+ totalEDof);
                                            
                        totalvDofBF_[j][k] = 0.0;
                        
                        Cv_pBoundary[j][k] = CvBoundary[j][k]/6.02214129e23;
                        
                        scalar gasConstant = 0.0;
                        scalar gamma = 0.0;
                        scalar speedOfSound = 0.0;
                        
                        if(molMassBoundary[j][k] > VSMALL)
                        {
                            gasConstant =   // R = k/m
                               physicoChemical::k.value()/molMassBoundary[j][k];
                        }
                        
                        if(CvBoundary[j][k] > VSMALL)
                        {
                            gamma = CpBoundary[j][k]/CvBoundary[j][k];
                        }
                        
                        if(gamma > VSMALL && gasConstant > VSMALL && 
                            translationalT_.boundaryFieldRef()[j][k] > VSMALL) // ***
                        {
                            speedOfSound = sqrt
                                    (
                                            gamma*gasConstant*                   
                                        translationalT_.boundaryField()[j][k]
                                    );
                        }
                        
                        if(speedOfSound > VSMALL)
                        {
                            Ma_.boundaryFieldRef()[j][k] =
                                        mag(UMean_.boundaryFieldRef()[j][k])/
                                        speedOfSound; // ***
                        }
                        else
                        {
                            Ma_.boundaryFieldRef()[j][k] = 0.0; // ***
                        }
                        
                        q_.boundaryFieldRef()[j][k] = qBF_[j][k]/nAvTimeSteps; // ***
                        
                        fD_.boundaryFieldRef()[j][k] = fDBF_[j][k]/nAvTimeSteps; // ***
                    }

                    p_.boundaryFieldRef()[j] = fD_.boundaryField()[j] & n_[j]; // ***
                    
                    tau_.boundaryFieldRef()[j] = sqrt(
                        sqr(fD_.boundaryField()[j] & t1_[j])
                        + sqr(fD_.boundaryField()[j] & t2_[j])); // ***
                }
            }

            Pout << "boundaryCells_" << endl;
            forAll(boundaryCells_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                
                labelList bCs = boundaryCells_[j];
                
                forAll(bCs, k)
                {
                    if(isA<polyPatch>(patch))
                    {
                        if(!isA<emptyPolyPatch>(patch))
                        {
                            if(!isA<cyclicPolyPatch>(patch))
                            {
                                dsmcRhoNMean_.boundaryFieldRef()[j][k] =
                                                    dsmcRhoNMean_[bCs[k]];  // changed .boundaryField() to .boundaryFieldRef()
                                rhoN_.boundaryFieldRef()[j][k] =
                                                    rhoN_[bCs[k]];  // ***
                                rhoM_.boundaryFieldRef()[j][k] =
                                                    rhoM_[bCs[k]];  // ***
                                
                                if(measureMeanFreePath_)
                                {
                                    MFP_.boundaryFieldRef()[j][k] =
                                                        MFP_[bCs[k]];   // ***
                                    SOF_.boundaryFieldRef()[j][k] =
                                                        SOF_[bCs[k]];   // ***
//                                    MFPdX_.boundaryFieldRef()[j][k] =
//                                                        MFPdX_[bCs[k]]; // ***
                                    dXMFP_.boundaryFieldRef()[j][k] =
                                            dXMFP_[bCs[k]]; // ***
                                    MCR_.boundaryFieldRef()[j][k] =
                                                        MCR_[bCs[k]];   // ***
                                    MCT_.boundaryFieldRef()[j][k] =
                                                        MCT_[bCs[k]];   // ***
                                     MCTdt_.boundaryFieldRef()[j][k] =
                                                        MCTdt_[bCs[k]]; // ***
                                }
                                if(measureClassifications_)
                                {
                                    cIDist_.boundaryFieldRef()[j][k] =
                                                        cIDist_[bCs[k]];    // ***
                                    cIIDist_.boundaryFieldRef()[j][k] =
                                                        cIIDist_[bCs[k]];   // ***
                                    cIIIDist_.boundaryFieldRef()[j][k] =
                                                        cIIIDist_[bCs[k]];  // ***
                                }
                                if(measureHeatFluxShearStress_)
                                {                                          
                                    shearStressTensor_.boundaryFieldRef()[j][k] =
                                                    shearStressTensor_[bCs[k]]; // ***
                                    heatFluxVector_.boundaryFieldRef()[j][k] =
                                                        heatFluxVector_[bCs[k]];    // ***
                                    pressureTensor_.boundaryFieldRef()[j][k] =
                                                        pressureTensor_[bCs[k]];    // ***
                                }
                                            
                                if(!isA<wallPolyPatch>(patch))
                                {    
                                    translationalT_.boundaryFieldRef()[j][k] =
                                                    translationalT_[bCs[k]];    // ***
                                    rotationalT_.boundaryFieldRef()[j][k] =
                                                    rotationalT_[bCs[k]];   // ***
                                    vibrationalT_.boundaryFieldRef()[j][k] =
                                                    vibrationalT_[bCs[k]];  // ***
                                    overallT_.boundaryFieldRef()[j][k] =
                                                    overallT_[bCs[k]];  // ***
                                    p_.boundaryFieldRef()[j][k] =
                                                    p_[bCs[k]]; // ***
                                    Ma_.boundaryFieldRef()[j][k] =
                                                    Ma_[bCs[k]];    // ***
                                    UMean_.boundaryFieldRef()[j][k] =
                                                    UMean_[bCs[k]]; // ***
                                }
                                
                                
                            }
                        }
                    }
                }
//                 if(measureMeanFreePath_)
//                 {
//                     if(!isA<emptyPolyPatch>(patch))
//                     {
//                         if(!isA<cyclicPolyPatch>(patch))
//                         {
//                             forAll(boundaryCells_[j], k)
//                             {
//                                 MFP_.boundaryField()[j][k] = 
//MFP_[boundaryCells_[j][k]];
//                                 SOF_.boundaryField()[j][k] = 
//SOF_[boundaryCells_[j][k]];
//                                 MFPdX_.boundaryField()[j][k] = 
//MFPdX_[boundaryCells_[j][k]];
//                                 MCR_.boundaryField()[j][k] = 
//MCR_[boundaryCells_[j][k]];
//                                 MCT_.boundaryField()[j][k] = 
//MCT_[boundaryCells_[j][k]];
//                                 
//MCTdt_.boundaryField()[j][k] = 
//MCTdt_[boundaryCells_[j][k]];
//                             }
//                         }
//                     }
//                 }
//                 if(measureClassifications_)
//                 {
//                     if(isA<polyPatch>(patch))
//                     {
//                         if(!isA<emptyPolyPatch>(patch))
//                         {
//                             if(!isA<cyclicPolyPatch>(patch))
//                             {
//                                 forAll(boundaryCells_[j], k)
//                                 {
//                                     cIDist_.boundaryField()[j][k] = 
// cIDist_[boundaryCells_[j][k]];
//                                     cIIDist_.boundaryField()[j][k] = 
// cIIDist_[boundaryCells_[j][k]];
//                                     cIIIDist_.boundaryField()[j][k] 
// = cIIIDist_[boundaryCells_[j][k]];
//                                 }
//                             }
//                         }
//                     }
//                 }
//                 if(measureHeatFluxShearStress_)
//                 { 
//                     if(isA<polyPatch>(patch))
//                     {
//                         if(!isA<emptyPolyPatch>(patch))
//                         {
//                             if(!isA<cyclicPolyPatch>(patch))
//                             {
//                                 forAll(boundaryCells_[j], k)
//                                 {
//                                     shearStressTensor_.boundaryField()[j][k] 
//= shearStressTensor_[boundaryCells_[j][k]];
//                                     heatFluxVector_.boundaryField()[j][k] = 
// heatFluxVector_[boundaryCells_[j][k]];
//                                     pressureTensor_.boundaryField()[j][k] = 
// pressureTensor_[boundaryCells_[j][k]];
//                                 }
//                             }
//                         }
//                     }
//                 }
            }

            /*
            Pout << "  write out now" << endl;
            if(measureMeanFreePath_)
            {
                MFP_.write();
                //MFPdX_.write();
                dXMFP_.write();
                MCR_.write();
                MCT_.write();
                MCTdt_.write();
                MCS_.write();
                SOF_.write();
            }
            
            if(measureClassifications_)
            {
                cIDist_.write();
                cIIDist_.write();
                cIIIDist_.write();
            }
            
            if(measureErrors_)
            {
                densityError_.write();
                velocityError_.write();
                temperatureError_.write();
                pressureError_.write();
            }
            
            if(measureHeatFluxShearStress_)
            {
                heatFluxVector_.write();
                pressureTensor_.write();
                shearStressTensor_.write();
            }
            
            p_.write();
            translationalT_.write();
            rotationalT_.write();
            vibrationalT_.write();
            electronicT_.write();
            overallT_.write();
            q_.write();
            tau_.write();
            Ma_.write();
            UMean_.write();
            fD_.write();
            */
        }

        //- reset
        if(time_.resetFieldsAtOutput())
        {
            Pout << "Reset field after output!" << nl << endl;
            nTimeSteps_ = 0.0;
            
            forAll(rhoNMean_, c)
            {
                rhoNMean_[c] = scalar(0.0);
                rhoMMean_[c] = scalar(0.0);
                linearKEMean_[c] = scalar(0.0);
                momentumMean_[c] = vector::zero;
                rotationalEMean_[c] = scalar(0.0);
                rotationalDofMean_[c] = scalar(0.0);
                rhoNMeanInt_[c] = scalar(0.0);
                molsElec_[c] = scalar(0.0),
                nClassI_[c] = scalar(0.0);
                nClassII_[c] = scalar(0.0);
                nClassIII_[c] = scalar(0.0);
                collisionSeparation_[c] = scalar(0.0);
                nColls_[c] = scalar(0.0);
                muu_[c] = scalar(0.0);
                muv_[c] = scalar(0.0);
                muw_[c] = scalar(0.0);
                mvv_[c] = scalar(0.0);
                mvw_[c] = scalar(0.0);
                mww_[c] = scalar(0.0);
                mcc_[c] = scalar(0.0);
                mccu_[c] = scalar(0.0);
                mccv_[c] = scalar(0.0);
                mccw_[c] = scalar(0.0);
                eu_[c] = scalar(0.0);
                ev_[c] = scalar(0.0);
                ew_[c] = scalar(0.0);
                e_[c] = scalar(0.0);
                rhoNMeanXnParticle_[c] = scalar(0.0);
                rhoMMeanXnParticle_[c] = scalar(0.0);
                momentumMeanXnParticle_[c] = vector::zero;
                linearKEMeanXnParticle_[c] = scalar(0.0);
            }
            
            forAll(electronicETotal_, iD)
            {
                forAll(electronicETotal_[iD], cell)
                {
                    electronicETotal_[iD][cell] = 0.0;
                    mccSpecies_[iD][cell] = 0.0;
                    nParcels_[iD][cell] = 0.0;
                    nGroundElectronicLevel_[iD][cell] = 0.0;
                    nFirstElectronicLevel_[iD][cell] = 0.0;
                    nParcelsXnParticle_[iD][cell] = 0.0;
                    
                    forAll(vibrationalETotal_[iD], v)
                    {
                       vibrationalETotal_[iD][v][cell] = 0.0; 
                    }
                }
            }
            
            // reset boundary information
            
            forAll(rhoNBF_, j)
            {
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                linearKEBF_[j] = 0.0;
                speciesRhoNIntBF_[j] = 0.0;
                speciesRhoNElecBF_[j] = 0.0;
                rotationalEBF_[j] = 0.0;
                rotationalDofBF_[j] = 0.0;
                qBF_[j] = 0.0;
                fDBF_[j] = vector::zero;
                momentumBF_[j] = vector::zero;
            }
            
            forAll(speciesRhoNBF_, i)
            {
                forAll(speciesRhoNBF_[i], j)
                { 
                    speciesRhoNBF_[i][j] = 0.0;
                    vibrationalEBF_[i][j] = 0.0;
                    electronicEBF_[i][j] = 0.0;
                    mccSpeciesBF_[i][j] = 0.0;
                }
            }
            Pout << "Reset done." << endl;
        }
        
        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }
        
    }
}


//- write field
void dsmcVolFields::writeField()
{
}

void dsmcVolFields::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}




} // End namespace Foam

// ************************************************************************* //

