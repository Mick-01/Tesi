/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "frozenSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void frozenSST<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut
    (
        S2
    );

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void frozenSST<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
frozenSST<BasicTurbulenceModel>::frozenSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    R_
    (
        IOobject
        (
            "R",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", sqr(dimLength)/sqr(dimTime)/dimTime, 0.0)
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

template<class BasicTurbulenceModel>
void frozenSST<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    volScalarField& R = this->R_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicTurbulenceModel::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(this->S2(tgradU()));
    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);


    // - boundary condition changes a cell value
    // - normally this would be triggered through correctBoundaryConditions
    // - which would do
    //      - fvPatchField::evaluate() which calls
    //      - fvPatchField::updateCoeffs()
    // - however any processor boundary conditions already start sending
    //   at initEvaluate so would send over the old value.
    // - avoid this by explicitly calling updateCoeffs early and then
    //   only doing the boundary conditions that rely on initEvaluate
    //   (currently only coupled ones)

    //- 1. Explicitly swap values on coupled boundary conditions
    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();
    // omegaWallFunctions change the cell value! Make sure to push these to
    // coupled neighbours. Note that we want to avoid the re-updateCoeffs
    // of the wallFunctions so make sure to bypass the evaluate on
    // those patches and only do the coupled ones.
    this->omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    ////- 2. Make sure the boundary condition calls updateCoeffs from
    ////     initEvaluate
    ////     (so before any swap is done - requires all coupled bcs to be
    ////      after wall bcs. Unfortunately this conflicts with cyclicACMI)
    //this->omega_.correctBoundaryConditions();


    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    const volScalarField F1(this->F1(CDkOmega));
    const volScalarField F23(this->F23());

    {
        // Turbulent kinetic energy equation
        
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, this->k_)
          + fvm::div(alphaRhoPhi, this->k_)
          - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
         ==
            alpha()*rho()*this->Pk(G)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
          - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, tgradU()), this->k_)
          + alpha()*rho()*this->betaStar_*this->omegaInf_*this->kInf_
          + this->kSource()
          + fvOptions(alpha, rho, this->k_)
        );

        tgradU.clear();

        kEqn.ref().relax();
        fvOptions.constrain(kEqn.ref());
        //solve(kEqn);
        R.primitiveFieldRef() = kEqn().residual();
        fvOptions.correct(this->k_);
        bound(this->k_, this->kMin_);

        /*R = fvc::ddt(alpha, rho, this->k_)
        + fvc::div(alphaRhoPhi, this->k_)
        - fvc::laplacian(alpha*rho*this->DkEff(F1), this->k_)
        - (
            alpha()*rho()*this->Pk(G)
          - fvc::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
          - fvc::Sp(alpha()*rho()*this->epsilonByk(F1, tgradU()), this->k_)
          + alpha()*rho()*this->betaStar_*this->omegaInf_*this->kInf_
          + this->kSource()
          + fvOptions(alpha, rho, this->k_)
        );*/
    }

    {
        const volScalarField::Internal gamma(this->gamma(F1));
        const volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = this->GbyNu(GbyNu0, F23(), S2());

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          + fvm::div(alphaRhoPhi, this->omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
            alpha()*rho()*gamma*(GbyNu0+R.internalField()/(nut))
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
          - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
          + alpha()*rho()*beta*sqr(this->omegaInf_)
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvOptions(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    correctNut(S2);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
