
tmp<volScalarField> blend
(
    const volScalarField& F1,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
) const
{
    return F1*(psi1 - psi2) + psi2;
}
tmp<volScalarField::Internal> blend
(
    const volScalarField& F1,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
) const
{
    return F1*(psi1 - psi2) + psi2;
}

tmp<volScalarField> F1 () const{
    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k)/(omega*y),
                scalar(500)*(nu)/(sqr(y)*omega)
            ),
            (4*sigmaOmega2)*k/(CDkOmega*sqr(y))
        ),
        scalar(10)
    );
    return tanh(pow4(arg1));
}
// tmp<volScalarField> alphaK(const volScalarField& F1) const
// {
//     return blend(F1, alphaK1_, alphaK2_);
// }

// tmp<volScalarField> alphaOmega(const volScalarField& F1) const
// {
//     return blend(F1, alphaOmega1_, alphaOmega2_);
// }

// tmp<volScalarField::Internal> beta
// (
//     const volScalarField::Internal& F1
// ) const
// {
//     return tmp<volScalarField::Internal>::New
//     (
//         IOobject::scopedName(this->type(), "beta"),
//         blend(F1, beta1_, beta2_)
//     );
// }

// tmp<volScalarField::Internal> gamma
// (
//     const volScalarField::Internal& F1
// ) const
// {
//     return tmp<volScalarField::Internal>::New
//     (
//         IOobject::scopedName(this->type(), "gamma"),
//         blend(F1, gamma1_, gamma2_)
//     );
// }
