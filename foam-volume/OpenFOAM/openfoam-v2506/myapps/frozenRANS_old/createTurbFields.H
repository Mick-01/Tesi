IOdictionary transportDict
(
    IOobject
    (
        "TransportProperties",      
        runTime.constant(),  
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar& nu_
(
    transportDict.lookup("nu")
);

volScalarField omega_ 
(
    IOobject
    (
        "omega",                  
        runTime.timeName(),         
        mesh,                       
        IOobject::READ_IF_PRESENT, 
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("omega", dimless/dimTime,1.0)
);

volScalarField nut_ 
(
    IOobject
    (
        "nut",                  
        runTime.timeName(),         
        mesh,                       
        IOobject::NO_READ, 
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("nut", sqr(dimLength)/dimTime,1.0)
);

volScalarField R_
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    dimensionedScalar("zero", k.dimensions()/dimTime, 0.0)
);
