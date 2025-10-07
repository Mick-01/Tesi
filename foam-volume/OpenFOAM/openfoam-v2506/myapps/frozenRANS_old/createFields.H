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
/*volScalarField Omega
(
    IOobject
    (
        "omega",                  // field name
        runTime.timeName(),         // time directory
        mesh,                       // mesh database
        IOobject::MUST_READ_IF_MODIFIED,          //::NO_READ don’t read from disk
        IOobject::NO_WRITE          // don’t write to disk
    ),
    mesh,
    // dimensionedScalar("omega", dimless/dimTime, 1.0)
);
*/
