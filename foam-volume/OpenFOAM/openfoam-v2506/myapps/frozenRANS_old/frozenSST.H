#ifndef FROZENSST_H
#define FROZENSST_H
#include "wallDist.H"
#include "fvCFD.H"
#include "bound.H"
using namespace Foam;
class frozenSST {
    private:
        
        // Passed from external call
        const volVectorField& U_;
        const volScalarField& k_;
        const fvMesh& mesh_;

        // read from file
        dimensionedScalar nu_;
        

        // instant calculation
        volScalarField y_;
        volSymmTensorField S_;
        surfaceScalarField phi_;
        volScalarField divU_;

        // standard parameters k-omega SST
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

        // blended parameters, updated once per frame
        volScalarField sigmaK_;
        volScalarField sigmaOmega_;
        volScalarField gamma_;
        volScalarField beta_;

        // target, kept to pass information between times
        volScalarField omega_;
        volScalarField nut_;
        volScalarField R_;
    protected:
        // generates the blended parameters
        // @param F1 the field result of F1 function, to update once per cycle
        // @param psi1 the k-omega parameter
        // @param psi2 the k-epsilon parameter
        // @return the blended parameter field
        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        );

        // updates all the blended parameters 
        void updateBlended();

        // generates the difference term between k-omega and k-epsilon 
        tmp<volScalarField> CDKOmega();

        // returns the scalar field of the blending function F1
        // @params CDKOmega_ the CDKOmega value for the iteration
        tmp<volScalarField> F1 (const volScalarField& CDKOmega_);
        
        // returns a blending function for the nu_t equation
        tmp<volScalarField> F2 ();

        // returns the turbulent viscosity field
        tmp<volScalarField> nut ();

    public:
        // constructor, needs only U, k and the mesh;
        // the other parameters that need to be passed will be read directly
        frozenSST(
            const volVectorField& U,
            const volScalarField& k,
            const fvMesh& mesh
        );

        // a small utility to debug the class, variable usage
        void debug();

        // performs a single iteration step
        void cycle();
};
#endif
