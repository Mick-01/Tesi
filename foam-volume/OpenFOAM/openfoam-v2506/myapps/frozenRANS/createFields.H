volVectorField U (
    IOobject
    (
        "U",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    ),
    mesh
);

#include "createPhi.H"

singlePhaseTransportModel laminarTransport(U, phi);

volScalarField k (
    IOobject
    (
        "k",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField y(
    "y",
    wallDist::New(mesh).y()
);

IOdictionary transportProperties(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu (
    "nu",
    dimArea/dimTime,
    transportProperties
);
// SST params
dimensionedScalar sigmaK1       (0.85);
dimensionedScalar sigmaK2       (1.0);

dimensionedScalar sigmaOmega1   (0.5);
dimensionedScalar sigmaOmega2   (0.856);

dimensionedScalar gamma1        (5/9);
dimensionedScalar gamma2        (0.44);

dimensionedScalar beta1         (0.075);
dimensionedScalar beta2         (0.0828);

dimensionedScalar betaStar      (0.09);

dimensionedScalar a1            (0.31);
dimensionedScalar b1            (1.0);
dimensionedScalar c1            (10.0);


// Target fields

volScalarField omega (
    IOobject
    (
        "omega",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField nut (
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("nut", sqr(dimLength)/dimTime, 1.0)
);

// Helper fields

volScalarField CDkOmega(
    IOobject
    (
        "CDkOmega",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cdkomega", dimless/sqr(dimTime), 1.0)
);

volScalarField F1(
    IOobject
    (
        "F1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(1.0)
);

volScalarField sigmaK(
    IOobject
    (
        "sigmaK",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    sigmaK1
);

volScalarField sigmaOmega(
    IOobject
    (
        "sigmaomega",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    sigmaOmega1
);

volScalarField gamma(
    IOobject
    (
        "gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    gamma1
);

volScalarField beta(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    beta1
);

volScalarField F2(
    IOobject
    (
        "F2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(1.0)
);

volScalarField Pk(
    IOobject
    (
        "Pk",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Pk",sqr(dimLength)/sqr(dimTime)/dimTime,1.0)
);

volScalarField R(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R",sqr(dimLength)/sqr(dimTime)/dimTime,1.0)
);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
