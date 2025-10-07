#ifndef FROZENSST_H
#define FROZENSST_H
#include "wallDist.H"
#include "fvCFD.H"
using namespace Foam;
class frozenSST {
    private:
        volScalarField omega_;  //Target

        const volScalarField& k_;
        const volVectorField& U_;

        const dimensionedScalar& nu_;
        const fvMesh& mesh_;

        const volScalarField& y_;

        tmp<volScalarField> sigmaK_;
        tmp<volScalarField> sigmaOmega_;
        tmp<volScalarField> gamma_;
        tmp<volScalarField> beta_;

        dimensionedScalar sigmaK1_  = 0.85;
        dimensionedScalar sigmaK2_  = 1.0;

        dimensionedScalar sigmaOmega1_ = 0.5;
        dimensionedScalar sigmaOmega2_ = 0.856;

        dimensionedScalar gamma1_ = 5/9;
        dimensionedScalar gamma2_ = 0.44;

        dimensionedScalar beta1_ = 0.075;
        dimensionedScalar beta2_  = 0.0828;

        dimensionedScalar betaStar_ = 0.09;

        dimensionedScalar a1_ = 0.31;
        dimensionedScalar b1_ = 1.0;
        dimensionedScalar c1_ = 10.0;
    protected:
        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        );

        void updateBlended();
        tmp<volScalarField> CDKOmega();

    public:
        frozenSST(
            volScalarField& omega,
            const volScalarField& k,
            const volVectorField& U,
            const dimensionedScalar& nu,
            const fvMesh& mesh
        );

        tmp<volScalarField> F1 (const volScalarField& CDKOmega_);
        tmp<volScalarField> F2 ();
        tmp<volScalarField> nut_ ();
        void debug();
        void cycle();
};
#endif
