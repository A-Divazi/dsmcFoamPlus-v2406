/*---------------------------------------------------------------------------*\
  User function object: DSMC statistical uncertainty (implementation).
\*---------------------------------------------------------------------------*/

#include "dsmcStatistics.H"
#include "volFields.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace Foam::constant;

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(dsmcStatistics, 0);
    addToRunTimeSelectionTable(functionObject, dsmcStatistics, dictionary);
}
}


namespace Foam
{

static scalar tCrit95(const scalar df)
{
    if (df >= 30.0) return 1.96;
    static const scalar tbl[] =
    {
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
        2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
        2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042
    };
    const label i = clamp(label(df), label(1), label(30)) - 1;
    return tbl[i];
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::functionObjects::dsmcStatistics::dsmcStatistics
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    equilibrationTime_(0.0),
    sampleInterval_(10),
    reportInterval_(100),
    confidenceLevel_(0.95),
    uncertaintyTolerance_(0.02),
    driftTolerance_(0.01),
    minimumEffectiveSamples_(100),
    minimumSamples_(20),
    blockCount_(10),
    maxLag_(200),
    requiredConsecutivePasses_(3),
    resetOnRestart_(true),
    step_(0),
    equilSteps_(0),
    samples_(0),
    series_(6),
    prevUrel_(6, -1.0),
    latchedPass_(6, false),
    consecPasses_(0),
    converged_(false),
    datasetPrefix_("postProcessing/statisticalUncertainty")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::dsmcStatistics::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("equilibrationTime", equilibrationTime_);
    dict.readIfPresent("sampleInterval", sampleInterval_);
    dict.readIfPresent("reportInterval", reportInterval_);
    dict.readIfPresent("confidenceLevel", confidenceLevel_);
    dict.readIfPresent("uncertaintyTolerance", uncertaintyTolerance_);
    dict.readIfPresent("driftTolerance", driftTolerance_);
    dict.readIfPresent("minimumEffectiveSamples", minimumEffectiveSamples_);
    dict.readIfPresent("minimumSamples", minimumSamples_);
    dict.readIfPresent("blockCount", blockCount_);
    dict.readIfPresent("maxLag", maxLag_);
    dict.readIfPresent("requiredConsecutivePasses", requiredConsecutivePasses_);
    dict.readIfPresent("resetOnRestart", resetOnRestart_);

    sampleInterval_ = max(sampleInterval_, label(1));
    reportInterval_ = max(reportInterval_, sampleInterval_);

    return true;
}


void Foam::functionObjects::dsmcStatistics::sample()
{
    const auto* rhoMP = mesh_.cfindObject<volScalarField>("rhoM");
    const auto* momP = mesh_.cfindObject<volVectorField>("momentum");
    const auto* lkP = mesh_.cfindObject<volScalarField>("linearKE");
    const auto* rhoNP = mesh_.cfindObject<volScalarField>("rhoN");

    if (!rhoMP || !momP || !lkP || !rhoNP)
    {
        return;
    }

    const volScalarField& rhoM = *rhoMP;
    const volVectorField& momentum = *momP;
    const volScalarField& linKE = *lkP;
    const volScalarField& rhoN = *rhoNP;

    const scalarField& V = mesh_.V();

    scalar rhoS = 0, TAv = 0, pAv = 0, volS = 0, massS = 0;
    vector momS = Zero;

    forAll(V, celli)
    {
        const scalar r = max(rhoM[celli], SMALL);
        const vector U = momentum[celli]/r;
        const scalar lk = max(linKE[celli], scalar(0));
        const scalar n = max(rhoN[celli], SMALL);
        const scalar trans = max(lk - 0.5*r*(U & U), scalar(0));
        const scalar Tc = 2.0/(3.0*physicoChemical::k.value()*n)*trans;
        const scalar pc = physicoChemical::k.value()*n*Tc;
        const scalar v = V[celli];
        rhoS += v*r;
        momS += v*momentum[celli];
        massS += v*r;
        TAv += v*Tc;
        pAv += v*pc;
        volS += v;
    }

    reduce(rhoS, sumOp<scalar>());
    reduce(momS, sumOp<vector>());
    reduce(massS, sumOp<scalar>());
    reduce(TAv, sumOp<scalar>());
    reduce(pAv, sumOp<scalar>());
    reduce(volS, sumOp<scalar>());

    const scalar rho = (volS > 0) ? rhoS/volS : 0;
    const vector U = (massS > 0) ? momS/massS : Zero;
    const scalar T = (volS > 0) ? TAv/volS : 0;
    const scalar p = (volS > 0) ? pAv/volS : 0;

    series_[0].append(rho);
    series_[1].append(U.x());
    series_[2].append(U.y());
    series_[3].append(U.z());
    series_[4].append(T);
    series_[5].append(p);
    ++samples_;

    if (Pstream::master())
    {
        const fileName dir(mesh_.time().globalPath()/datasetPrefix_);
        mkDir(dir);
        std::ofstream os((dir/"raw_samples.csv").c_str(), std::ios::app);
        os << mesh_.time().value() << ' ' << rho << ' '
           << U.x() << ' ' << U.y() << ' ' << U.z() << ' '
           << T << ' ' << p << std::endl;
    }
}


// Compute one quantity's statistics from its time series.
static void computeStats
(
    const scalarList& x,
    const scalar scale,
    const scalar driftTol,
    const scalar uncTol,
    const label minEff,
    const label minSamples,
    const label maxLag,
    const scalar prevUrel,
    scalar& mean, scalar& sd, scalar& tauInt, scalar& Neff, scalar& SE,
    scalar& Urel, scalar& drift, scalar& CIlo, scalar& CIhi, bool& pass
)
{
    const label N = x.size();
    mean = sd = tauInt = Neff = SE = Urel = drift = 0;
    CIlo = CIhi = 0;
    pass = false;
    if (N < 1)
    {
        return;
    }

    for (const scalar v : x) mean += v;
    mean /= N;

    if (N > 1)
    {
        scalar s2 = 0;
        for (const scalar v : x) s2 += sqr(v - mean);
        sd = std::sqrt(s2/(N - 1));
    }

    const label maxK = min(maxLag, N - 1);
    scalar C0 = 0;
    for (const scalar v : x) C0 += sqr(v - mean);
    C0 /= N;

    tauInt = 0.5;
    if (C0 > 0)
    {
        for (label k = 1; k <= maxK; ++k)
        {
            scalar Ck = 0;
            for (label i = 0; i < N - k; ++i) Ck += (x[i] - mean)*(x[i + k] - mean);
            Ck /= (N - k);
            const scalar rk = Ck/C0;
            if (rk <= 0.0 || rk < 1e-6) break;
            tauInt += rk;
        }
    }

    Neff = (tauInt > 0) ? N/(2.0*tauInt) : N;
    SE = (Neff >= 1) ? sd/std::sqrt(max(Neff, 1.0)) : 0;

    const scalar eps = 1e-12;
    const scalar scaleMag = max(mag(scale), eps);
    Urel = 100.0*SE/scaleMag;

    if (N >= 8)
    {
        const label h = N/2;
        scalar m1 = 0, m2 = 0;
        for (label i = 0; i < h; ++i) m1 += x[i];
        for (label i = h; i < N; ++i) m2 += x[i];
        m1 /= h; m2 /= (N - h);
        drift = mag(m2 - m1)/scaleMag;
    }

    const scalar t = tCrit95(max(Neff - 1.0, 1.0));
    CIlo = mean - t*SE;
    CIhi = mean + t*SE;

    const bool effSamples = (N >= minSamples);
    const bool A = drift < driftTol;
    const bool B = Urel < 100.0*uncTol;
    const bool C = (Neff >= minEff);
    const bool D = (prevUrel < 0) ||
                   (mag(Urel - prevUrel) < 0.10*max(mag(prevUrel), 1e-9));
    pass = (effSamples && A && B && C && D);
}

} // End namespace Foam


bool Foam::functionObjects::dsmcStatistics::execute()
{
    ++step_;
    const scalar dt = mesh_.time().deltaTValue();
    equilSteps_ = (dt > 0) ? label(equilibrationTime_/dt + 0.5) : 0;

    const label rel = step_ - equilSteps_;

    if (rel > 0 && (rel % sampleInterval_ == 0))
    {
        sample();
    }

    if (rel > 0 && (rel % reportInterval_ == 0))
    {
        writeStats();
        writeSummary();
    }

    return true;
}


bool Foam::functionObjects::dsmcStatistics::write()
{
    return true;
}


void Foam::functionObjects::dsmcStatistics::writeStats()
{
    static const wordList names = { "rho", "Ux", "Uy", "Uz", "T", "p" };

    bool allPass = true;
    // Ux mean used as the characteristic scale for the near-zero transverse
    // velocities (Uy, Uz) so their relative uncertainty/drift is meaningful.
    scalar meanUx = 0;
    {
        const scalarList& ux = series_[1];
        for (const scalar v : ux) meanUx += v;
        meanUx = (ux.size()) ? meanUx/ux.size() : 0;
    }
    for (label q = 0; q < names.size(); ++q)
    {
        scalarList x(series_[q].size());
        forAll(x, i) x[i] = series_[q][i];

        scalar scale = meanUx;
        if (q != 2 && q != 3)
        {
            scale = 0;
            for (const scalar v : x) scale += v;
            scale = (x.size()) ? scale/x.size() : 0;
        }

        scalar mean, sd, tauInt, Neff, SE, Urel, drift, CIlo, CIhi;
        bool pass = false;
        computeStats(x, scale, driftTolerance_, uncertaintyTolerance_,
                     minimumEffectiveSamples_, minimumSamples_, maxLag_,
                     prevUrel_[q], mean, sd, tauInt, Neff, SE, Urel, drift,
                     CIlo, CIhi, pass);

        prevUrel_[q] = Urel;
        // latch: once a quantity passes it stays PASS
        if (pass)
        {
            latchedPass_[q] = true;
        }
        const bool qPass = latchedPass_[q];
        allPass = allPass && qPass;

        if (Pstream::master())
        {
            const fileName dir(mesh_.time().globalPath()/datasetPrefix_);
            mkDir(dir);
            const fileName csv(dir/"statistics.csv");
            const bool first = !isFile(csv);
            std::ofstream os(csv.c_str(), std::ios::app);
            if (first)
            {
                os << "#time quantity location N mean stdDev tauInt effectiveSamples "
                      "standardError relativeUncertainty CI95Lower CI95Upper drift converged\n";
            }
            os << mesh_.time().value() << ' ' << names[q] << " domain "
               << x.size() << ' ' << mean << ' ' << sd << ' ' << tauInt << ' '
               << Neff << ' ' << SE << ' ' << Urel << ' ' << CIlo << ' ' << CIhi
               << ' ' << drift << ' ' << (qPass ? 1 : 0) << std::endl;
        }
    }

    const bool thisPass = allPass;
    consecPasses_ = thisPass ? (consecPasses_ + 1) : 0;
    converged_ = (consecPasses_ >= requiredConsecutivePasses_);

    if (thisPass && Pstream::master())
    {
        Log << type() << ' ' << name() << " : all quantities pass (consecutive "
            << consecPasses_ << '/' << requiredConsecutivePasses_ << ')'
            << (converged_ ? " -> CONVERGED" : "") << nl;
    }
}


void Foam::functionObjects::dsmcStatistics::writeSummary()
{
    if (!Pstream::master())
    {
        return;
    }
    static const wordList names = { "rho", "Ux", "Uy", "Uz", "T", "p" };

    const fileName dir(mesh_.time().globalPath()/datasetPrefix_);
    mkDir(dir);
    std::ofstream os((dir/"summary.dat").c_str(), std::ios::trunc);

    os << "# t=" << mesh_.time().value()
       << "  samples=" << samples_
       << "  converged=" << (converged_ ? "yes" : "no")
       << "  consecutivePasses=" << consecPasses_
       << "/" << requiredConsecutivePasses_ << "\n";
    os << std::left << std::setw(9) << "quantity"
       << std::right << std::setw(15) << "mean"
       << std::setw(15) << "SE"
       << std::setw(10) << "Urel%"
       << std::setw(9) << "Neff"
       << std::setw(9) << "tauInt"
       << std::setw(13) << "drift"
       << ' ' << std::left << std::setw(6) << "status" << "\n";
    scalar meanUx = 0;
    {
        const scalarList& ux = series_[1];
        for (const scalar v : ux) meanUx += v;
        meanUx = (ux.size()) ? meanUx/ux.size() : 0;
    }
    for (label q = 0; q < names.size(); ++q)
    {
        scalarList x(series_[q].size());
        forAll(x, i) x[i] = series_[q][i];
        scalar scale = meanUx;
        if (q != 2 && q != 3)
        {
            scale = 0;
            for (const scalar v : x) scale += v;
            scale = (x.size()) ? scale/x.size() : 0;
        }
        scalar mean, sd, tauInt, Neff, SE, Urel, drift, CIlo, CIhi;
        bool pass = false;
        computeStats(x, scale, driftTolerance_, uncertaintyTolerance_,
                     minimumEffectiveSamples_, minimumSamples_, maxLag_,
                     prevUrel_[q], mean, sd, tauInt, Neff, SE, Urel, drift,
                     CIlo, CIhi, pass);
        os << std::left << std::setw(9) << names[q]
           << std::right << std::setprecision(5) << std::scientific
           << std::setw(15) << mean << std::setw(15) << SE
           << std::setprecision(4) << std::fixed
           << std::setw(10) << Urel
           << std::setw(9) << Neff << std::setw(9) << tauInt
           << std::setprecision(3) << std::scientific << std::setw(13) << drift
           << ' ' << std::left << std::setw(6)
           << (latchedPass_[q] ? "PASS" : "FAIL") << "\n";
    }
}
