/*---------------------------------------------------------------------------*\
  User function object: DSMC statistical uncertainty / convergence monitoring.

  Samples the instantaneous domain-average of rho, Ux, Uy, Uz, T, p after an
  equilibration period, and each report interval computes:
    mean, std, autocorrelation -> tauInt -> N_eff = N/(2*tauInt) (NOT s/sqrt(N)),
    standard error SE = s/sqrt(N_eff), relative uncertainty, Student-t 95% CI,
    block-average uncertainty, and stationarity drift.

  Convergence requires (for EVERY quantity, no majority voting):
    A drift < driftTolerance
    B Urel < uncertaintyTolerance
    C N_eff >= minimumEffectiveSamples
    D |Urel change| < 10% between successive reports
  sustained for requiredConsecutivePasses consecutive reports.

  Outputs to postProcessing/statisticalUncertainty/:
    raw_samples.csv, statistics.csv, summary.dat

  Usage:
    dsmcStatistics1
    {
        type            dsmcStatistics;
        libs            (libUserDsmcFoam);
        equilibrationTime 0;        // seconds of transient to discard
        sampleInterval   10;        // steps between samples
        reportInterval   100;       // steps between statistics reports
        confidenceLevel  0.95;
        uncertaintyTolerance 0.02;
        driftTolerance   0.01;
        minimumEffectiveSamples 100;
        minimumSamples   20;
        blockCount       10;
        maxLag           200;
        requiredConsecutivePasses 3;
        resetOnRestart   true;
    }
\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_dsmcStatistics_H
#define Foam_functionObjects_dsmcStatistics_H

#include "fvMeshFunctionObject.H"
#include "DynamicList.H"

namespace Foam
{
namespace functionObjects
{

class dsmcStatistics
:
    public fvMeshFunctionObject
{
    // Config
    scalar equilibrationTime_;
    label sampleInterval_;
    label reportInterval_;
    scalar confidenceLevel_;
    scalar uncertaintyTolerance_;
    scalar driftTolerance_;
    label minimumEffectiveSamples_;
    label minimumSamples_;
    label blockCount_;
    label maxLag_;
    label requiredConsecutivePasses_;
    bool resetOnRestart_;

    // State
    label step_;
    label equilSteps_;
    label samples_;
    List<DynamicList<scalar>> series_;   // 6 quantities
    scalarList prevUrel_;
    List<bool> latchedPass_;             // once a quantity passes it stays PASS
    label consecPasses_;
    bool converged_;
    word datasetPrefix_;

    //- Domain-average sample of the 6 quantities; append to series_
    void sample();
    //- Write statistics.csv (append per quantity)
    void writeStats();
    //- Write summary.dat (overwrite)
    void writeSummary();

    dsmcStatistics(const dsmcStatistics&) = delete;
    void operator=(const dsmcStatistics&) = delete;

public:
    TypeName("dsmcStatistics");

    dsmcStatistics(const word& name, const Time& runTime, const dictionary& dict);
    virtual ~dsmcStatistics() = default;

    virtual bool read(const dictionary&);
    virtual bool execute();
    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
