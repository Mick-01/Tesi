#include "frozenSST.H"
#include "kLowReWallFunctionFvPatchScalarField.H"
#include <string>
// frozenSST::frozenSST
// (
//     const volVectorField& U,
//     const volScalarField& k,
//     const fvMesh& mesh
// ):
//     U_(U),
//     k_(k),
//     mesh_(mesh),

//     y_(wallDist::New(mesh).y()),
//     S_(symm(fvc::grad(U))),
//     phi_(fvc::flux(U))
//     {
//         divU_ = fvc::div(phi_);
//         #include "createTurbFields.H"
//         updateBlended();
//     }
bool dbg = 1;
void print_dbg(std::string s) {
    if (dbg){
        Info<<s<<endl;
    };
}
frozenSST::frozenSST
(
    const volVectorField& U,
    const volScalarField& k,
    const fvMesh& mesh
)
:
    U_(U),
    k_(k),
    mesh_(mesh),

    nu_
    (
        "nu",
        dimViscosity,
        IOdictionary
        (
            IOobject
            (
                "transportProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    ),

    y_(wallDist::New(mesh).y()),
    S_(symm(fvc::grad(U))),
    phi_(fvc::flux(U)),
    divU_
    (
        IOobject
        (
            "divU",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(phi_)
    ),

    sigmaK_
    (
        IOobject
        (
            "sigmaK",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    sigmaOmega_
    (
        IOobject
        (
            "sigmaOmega",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    gamma_
    (
        IOobject
        (
            "gamma",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    beta_
    (
        IOobject
        (
            "beta",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("omega", dimless/dimTime, 1.0)
    ),

    nut_
    (
        IOobject
        (
            "nut",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("nut", sqr(dimLength)/dimTime, 1.0)
    ),

    R_
    (
        IOobject
        (
            "R",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", k.dimensions()/dimTime, 0.0)
    )

{
    print_dbg("in the constructor");
    updateBlended();
}


tmp<volScalarField> frozenSST::blend
(
    const volScalarField& F1,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
)
{
    return F1*(psi1 - psi2) + psi2;
}

void frozenSST::updateBlended(){
    auto CDKOmega_ = CDKOmega();
    auto F1_ = F1(CDKOmega_);
    sigmaK_ = blend(F1_,sigmaK1_,sigmaK2_);
    sigmaOmega_ = blend(F1_,sigmaOmega1_,sigmaOmega2_);
    gamma_ = blend(F1_,gamma1_,gamma2_);
    beta_ = blend(F1_,beta1_,beta2_);
}

tmp<volScalarField> frozenSST::CDKOmega () {
    return max (
        2*sigmaOmega2_/omega_*(fvc::grad(k_) & fvc::grad(omega_)),
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );
}

tmp<volScalarField> frozenSST::F1(
    const volScalarField& CDKOmega_
){
    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(nu_)/(sqr(y_)*omega_)
            ),
            (4*sigmaOmega2_)*k_/(CDKOmega_*sqr(y_))
        ),
        scalar(10)
    );
    return tanh(pow4(arg1));
}

tmp<volScalarField> frozenSST::F2()
{
        tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(nu_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

tmp<volScalarField> frozenSST::nut()
{
    return a1_*k_/max(a1_*omega_,F2()*mag(S_));
}

void frozenSST::debug(){
    auto cd = CDKOmega();
    volScalarField cdField("cd", cd); 
    cdField.write();
    volScalarField("F1",F1(cdField)).write();
}

void frozenSST::cycle()
{
    print_dbg("Starting the cycle");
    print_dbg("Calculating CD_{k-omega}");
    auto cdko_ = CDKOmega();
    print_dbg("Calculating F1");
    auto F1_ = F1(cdko_);
    print_dbg("Updating blended fields");
    updateBlended();
    print_dbg("Calculating F2");
    auto F2_ = F2();
    print_dbg("Calculating nu_t");
    nut_ = nut(); 
    print_dbg("Calculating P_k");
    volScalarField Pk_ = min(2.0*nut_*magSqr(S_),10*betaStar_*omega_*k_);
    print_dbg("Starting the k equation");
    tmp<fvScalarMatrix> kEqn
    (
        fvm::div(phi_, k_)
        - fvm::laplacian(sigmaK_*nut()+nu_, k_)
        ==
        Pk_
        - fvm::SuSp((2.0/3.0)*divU_, k_)
        - fvm::Sp(betaStar_*omega_, k_)
        //TODO: add omegaInf and kInf (0 in the base case)
        // +betaStar_*omegaInf_*kInf_
        //+ kSource()                   //should be also 0 (?)
        //+ fvOptions(alpha, rho, k_)   // this too (?)
    );
     print_dbg("Finding R");
    R_.primitiveFieldRef() = kEqn().residual();
     print_dbg("Starting the omega equation");
        Info << "before prod\n";
    volScalarField prod =gamma_*(Pk_+R_)/nut_; 
    Info << "before diff\n";
    auto diff = fvm::laplacian(sigmaOmega_*nut_+nu_, omega_);
     Info << "before conv\n";
    auto conv = fvm::div(phi_, omega_);
    Info << "before suSp term\n";
    auto su = fvm::SuSp((2.0/3.0)*gamma_*divU_, omega_);
    Info << "before sp\n";
    auto sp = fvm::Sp(beta_*omega_,omega_);
    Info << "before CD term\n";
    auto cd = fvm::SuSp
        (
            (F1_ - scalar(1))*cdko_/omega_,
            omega_
        );
    Info << "about to assemble\n";
    tmp<fvScalarMatrix> omegaEqn = conv - diff == prod - su-sp-cd;
    Info << "after assemble\n";
    /*tmp<fvScalarMatrix> omegaEqn
    (
        fvm::div(phi_, omega_)
        - fvm::laplacian(sigmaOmega_*nut_+nu_, omega_)
        ==
        gamma_*(Pk_+R_)/nut_
        - fvm::SuSp((2.0/3.0)*gamma_*divU_, omega_)
        - fvm::Sp(beta_*omega_, omega_)
        - fvm::SuSp
        (
            (F1_ - scalar(1))*cdko_/omega_,
            omega_
        )
        //+ beta_*sqr(omegaInf_)
        //+ Qsas(2*magSqr(S_), gamma_, beta_) // should be 0
        //+ omegaSource()
        //+ fvOptions(alpha, rho, omega_)
    );*/
         print_dbg("relaxations");
        omegaEqn.ref().relax();
        //fvOptions.constrain(omegaEqn.ref());  //shouldn't be needed
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
         print_dbg("Solving for Omega");
        solve(omegaEqn);
        //fvOptions.correct(omega_);            //neither should this
        //omegaMin should be a parameter
         print_dbg("Bounding omega");
        bound(omega_,dimensionedScalar("omegaMin",dimless/dimTime,1e-10));
         print_dbg("Calculating nu_t");
        nut_ = nut();
}
