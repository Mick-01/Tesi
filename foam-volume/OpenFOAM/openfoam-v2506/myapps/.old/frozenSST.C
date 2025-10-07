#include "frozenSST.H"

tmp<volScalarField> frozenSST::blend
(
    const volScalarField& F1,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
)
{
    return F1*(psi1 - psi2) + psi2;
}
/*tmp<volScalarField::Internal> frozenSST::blend
(
    const volScalarField& F1,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
) const
{
    return F1*(psi1 - psi2) + psi2;
}*/
frozenSST::frozenSST
(
    volScalarField& omega,
    const volScalarField& k,
    const volVectorField& U,
    const dimensionedScalar& nu,
    const fvMesh& mesh
):
    omega_(omega),
    k_(k),
    U_(U),
    nu_(nu),
    mesh_(mesh),
    y_(wallDist::New(mesh).y())
    {
        updateBlended();
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

tmp<volScalarField> frozenSST::nut_()
{
    auto S = mag(symm(fvc::grad(U_)));
    return a1_*k_/max(a1_*omega_,F2()*S);
}

void frozenSST::debug(){
    auto cd = CDKOmega();
    volScalarField cdField("cd", cd); 
    cdField.write();
    volScalarField("F1",F1(cdField)).write();
}

void frozenSST::cycle()
{
    tmp<surfaceScalarField> phi_ (fvc::flux(U_));
    tmp<volScalarField> 
    /* tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha_, k_)
        + fvm::div(alphaRhoPhi, k_)
        - fvm::laplacian(alpha*rho*DkEff(F1), k_)
        ==
        alpha()*rho()*Pk(G)
        - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
        - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
        + alpha()*rho()*betaStar_*omegaInf_*kInf_
        + kSource()
        + fvOptions(alpha, rho, k_)
    ); */
}
